#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.io as pio


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def extract_rank(lineage: str, prefix: str) -> str:
    if not isinstance(lineage, str):
        return "Unknown"
    for part in lineage.split(";"):
        part = part.strip()
        if part.startswith(prefix):
            val = part[len(prefix):]
            return val if val else "Unknown"
    return "Unknown"


def find_attrition_csvs(root: Path) -> list[Path]:
    return sorted(root.rglob("*_attrition_lineage*.csv"))


def build_dataframe(csv_path: Path, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    required = {
        "Lineage",
        "Total Count",
        "Mapped Read Proportion",
        "Unmapped Read Proportion",
        "Average GC Content",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} missing columns: {', '.join(sorted(missing))}")

    df = df.copy()
    df["sample_id"] = sample_id
    df["source_csv"] = str(csv_path)

    df["Total Count"] = pd.to_numeric(df["Total Count"], errors="coerce").fillna(0).astype(int)
    df["Mapped Read Proportion"] = pd.to_numeric(df["Mapped Read Proportion"], errors="coerce").fillna(0.0)
    df["Average GC Content"] = pd.to_numeric(df["Average GC Content"], errors="coerce").fillna(0.0)

    df["mapped_prop_pct"] = df["Mapped Read Proportion"] * 100.0
    df["gc_pct"] = df["Average GC Content"] * 100.0

    # Dot size proxy
    df["mapped_reads_est"] = df["Total Count"]

    # Color groups
    df["phylum"] = df["Lineage"].apply(lambda x: extract_rank(x, "p__"))

    return df


def plot_for_one_sample(df: pd.DataFrame, title: str):
    df = df.copy()

    # Bubble scaling (sqrt → better dynamic range)
    df["size_for_plot"] = (df["mapped_reads_est"] ** 0.5).clip(lower=1)

    fig = px.scatter(
        df,
        x="gc_pct",
        y="mapped_prop_pct",
        color="phylum",
        size="size_for_plot",
        size_max=60,   # 👈 bigger bubbles
        opacity=0.7,
        hover_name="Lineage",
        hover_data={
            "Lineage": True,
            "phylum": True,
            "Total Count": True,
            "mapped_reads_est": True,
            "gc_pct": ":.2f",
            "mapped_prop_pct": ":.2f",
            "sample_id": True,
            "source_csv": True,
            "size_for_plot": False,
        },
        labels={
            "gc_pct": "GC Content (%)",
            "mapped_prop_pct": "Mapped Read Proportion (%)",
            "phylum": "Phylum",
        },
        title=title,
        template="plotly_white",
        custom_data=["mapped_reads_est"],
    )

    fig.update_layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    dragmode="zoom",
    margin=dict(l=60, r=40, t=70, b=60),

    legend_title_text="Phylum",
    legend=dict(
        font=dict(size=16),          # ← legend labels bigger
        title=dict(font=dict(size=18)),  # ← legend title bigger
        itemsizing="constant"
        ),
    )

    fig.update_xaxes(range=[0, 100], showline=True, linewidth=2, mirror=True)
    fig.update_yaxes(range=[0, 100], showline=True, linewidth=2, mirror=True)

    fig.add_shape(
        type="rect",
        x0=0, x1=100,
        y0=0, y1=100,
        xref="x", yref="y",
        line=dict(width=3),
        fillcolor="rgba(0,0,0,0)",
        layer="above",
    )

    return fig


def write_clickable_html_with_slider(fig, out_html: Path):
    """
    Write HTML with:
      - click panel (lineage summary)
      - slider to hide small dots by mapped_reads_est
      - modebar stripped of pan tools to prevent infinite wandering
    """
    base_html = pio.to_html(
    fig,
    include_plotlyjs="cdn",
    full_html=True,
    config={
        "displayModeBar": True,
        "scrollZoom": True,
        "modeBarButtonsToRemove": ["select2d", "lasso2d"],

        # 👇 THIS is the key addition
        "toImageButtonOptions": {
            "format": "svg",     # default export format
            "filename": "attrition_plot",
            "height": 900,
            "width": 1400,
            "scale": 1
        }
    },
)

    panel = r"""
<div style="max-width: 1150px; margin: 0 auto; font-family: Arial, sans-serif;">
  <hr/>
  <div style="display:flex; gap: 18px; align-items: flex-start; flex-wrap: wrap;">
    <div style="flex: 1; min-width: 320px;">
      <h3 style="margin: 0.2rem 0;">Clicked point</h3>
      <div id="clicked-info" style="padding: 0.75rem; border: 1px solid #ddd; border-radius: 10px; background: #ffffff;">
        Click a dot to display its lineage + stats here. (Legend items toggle phyla.)
      </div>
    </div>

    <div style="width: 360px; min-width: 280px;">
      <h3 style="margin: 0.2rem 0;">Filter</h3>
      <div style="padding: 0.75rem; border: 1px solid #ddd; border-radius: 10px; background: #ffffff;">

        <div style="margin-bottom: 0.5rem;">
          <b>Minimum mapped reads (estimated):</b>
          <span id="minreads-label" style="font-family: monospace;">0</span>
        </div>
        <input id="minreads" type="range" min="0" max="1000" value="0" step="1" style="width: 100%;">
        <div style="font-size: 0.9rem; color: #333; margin-top: 0.4rem;">
          Dots below this threshold fade out and become unclickable.
        </div>

        <hr style="margin: 0.8rem 0; border: none; border-top: 1px solid #eee;">

        <div style="margin-bottom: 0.4rem;">
          <b>GC x-axis cutoffs (%):</b>
        </div>

        <div style="display:flex; gap: 10px; align-items: center; flex-wrap: wrap;">
          <label style="font-size: 0.9rem; color:#333;">
            Min:
            <input id="gc-min" type="number" min="0" max="100" step="0.1"
                   value="0" style="width: 90px; margin-left: 6px;">
          </label>

          <label style="font-size: 0.9rem; color:#333;">
            Max:
            <input id="gc-max" type="number" min="0" max="100" step="0.1"
                   value="100" style="width: 90px; margin-left: 6px;">
          </label>

          <button id="gc-apply"
                  style="padding: 0.35rem 0.6rem; border: 1px solid #ccc; border-radius: 8px; background: #fafafa; cursor: pointer;">
            Apply
          </button>

          <button id="gc-reset"
                  style="padding: 0.35rem 0.6rem; border: 1px solid #ccc; border-radius: 8px; background: #fafafa; cursor: pointer;">
            Reset
          </button>
        </div>

        <div id="gc-msg" style="font-size: 0.85rem; color: #666; margin-top: 0.4rem;">
          Tip: enter e.g. 25 and 75, then Apply.
        </div>

      </div>
    </div>
  </div>
</div>

<script>
(function() {
  // Find plot div
  var plotDiv = document.getElementsByClassName('plotly-graph-div')[0];
  if (!plotDiv) return;

  // Helper to compute a good slider max from data
  function computeMaxMappedReads() {
    var maxv = 0;
    (plotDiv.data || []).forEach(function(tr) {
      if (!tr.customdata) return;
      for (var i=0; i<tr.customdata.length; i++) {
        var v = tr.customdata[i] ? tr.customdata[i][0] : 0;
        if (v > maxv) maxv = v;
      }
    });
    return maxv;
  }

  // Clamp view: set axes back to 0–100 when user double-clicks
  plotDiv.on('plotly_doubleclick', function() {
    Plotly.relayout(plotDiv, {'xaxis.range':[0,100], 'yaxis.range':[0,100]});
    // also keep the input boxes in sync
    var gmin = document.getElementById('gc-min');
    var gmax = document.getElementById('gc-max');
    if (gmin) gmin.value = 0;
    if (gmax) gmax.value = 100;
  });

  // Click panel (with wrapping for long lineage)
  plotDiv.on('plotly_click', function(data) {
    if (!data || !data.points || !data.points.length) return;
    var p = data.points[0];

    var lineage = (p.hovertext || "NA");
    var x = (p.x !== undefined) ? p.x.toFixed(2) : "NA";
    var y = (p.y !== undefined) ? p.y.toFixed(2) : "NA";
    var mappedReads = (p.customdata && p.customdata[0] !== undefined) ? p.customdata[0] : "NA";

    var html = ""
      + "<div style='line-height:1.35;'>"
      + "<b>Lineage:</b> "
      + "<span style='display:inline; white-space:normal; overflow-wrap:anywhere; word-break:break-word;'>"
      + lineage
      + "</span>"
      + "<br/><b>GC (%):</b> " + x
      + "<br/><b>Mapped proportion (%):</b> " + y
      + "<br/><b>Mapped reads (est.):</b> " + mappedReads
      + "</div>";

    var box = document.getElementById('clicked-info');
    if (box) box.innerHTML = html;
  });

  // Slider filtering: fade out points with mapped_reads_est < threshold
  function applyThreshold(thr) {
    var update = { 'marker.opacity': [] };

    (plotDiv.data || []).forEach(function(tr) {
      var opac = [];
      if (tr.customdata) {
        for (var i=0; i<tr.customdata.length; i++) {
          var v = tr.customdata[i] ? tr.customdata[i][0] : 0;
          opac.push((v >= thr) ? 0.75 : 0.0);
        }
      }
      update['marker.opacity'].push(opac);
    });

    var traceIdxs = (plotDiv.data || []).map(function(_, i){ return i; });
    Plotly.restyle(plotDiv, update, traceIdxs);
  }

  // Init slider bounds
  var maxv = computeMaxMappedReads();
  var slider = document.getElementById('minreads');
  var label = document.getElementById('minreads-label');
  if (slider && label) {
    slider.max = String(Math.max(10, maxv));
    slider.value = "0";
    label.textContent = "0";

    slider.addEventListener('input', function() {
      var thr = parseInt(slider.value || "0", 10);
      label.textContent = String(thr);
      applyThreshold(thr);
    });
  }

  // --- GC cutoff controls ---
  function clamp(n, lo, hi) { return Math.max(lo, Math.min(hi, n)); }

  function applyGcRange(xmin, xmax) {
    Plotly.relayout(plotDiv, {'xaxis.range':[xmin, xmax]});
  }

  var gcMin = document.getElementById('gc-min');
  var gcMax = document.getElementById('gc-max');
  var gcApply = document.getElementById('gc-apply');
  var gcReset = document.getElementById('gc-reset');
  var gcMsg = document.getElementById('gc-msg');

  function setMsg(txt, isErr) {
    if (!gcMsg) return;
    gcMsg.textContent = txt;
    gcMsg.style.color = isErr ? "#b00020" : "#666";
  }

  function readGcInputs() {
    var xmin = parseFloat(gcMin ? gcMin.value : "0");
    var xmax = parseFloat(gcMax ? gcMax.value : "100");
    if (Number.isNaN(xmin) || Number.isNaN(xmax)) {
      setMsg("Please enter numeric GC min/max.", true);
      return null;
    }
    xmin = clamp(xmin, 0, 100);
    xmax = clamp(xmax, 0, 100);
    if (xmin >= xmax) {
      setMsg("GC min must be less than GC max.", true);
      return null;
    }
    // keep inputs synced to clamped values
    if (gcMin) gcMin.value = xmin;
    if (gcMax) gcMax.value = xmax;
    return [xmin, xmax];
  }

  if (gcApply) {
    gcApply.addEventListener('click', function() {
      var vals = readGcInputs();
      if (!vals) return;
      applyGcRange(vals[0], vals[1]);
      setMsg("Applied GC window: " + vals[0] + "–" + vals[1], false);
    });
  }

  // apply on Enter in either box
  function onEnter(e) {
    if (e.key === "Enter") {
      var vals = readGcInputs();
      if (!vals) return;
      applyGcRange(vals[0], vals[1]);
      setMsg("Applied GC window: " + vals[0] + "–" + vals[1], false);
    }
  }
  if (gcMin) gcMin.addEventListener('keydown', onEnter);
  if (gcMax) gcMax.addEventListener('keydown', onEnter);

  if (gcReset) {
    gcReset.addEventListener('click', function() {
      if (gcMin) gcMin.value = 0;
      if (gcMax) gcMax.value = 100;
      applyGcRange(0, 100);
      setMsg("Reset GC window to 0–100.", false);
    });
  }

  // Start with no filter
  applyThreshold(0);
})();
</script>
"""

    out = base_html.replace("</body>", panel + "\n</body>")
    out_html.write_text(out, encoding="utf-8")

def main():
    ap = argparse.ArgumentParser(
        description="Create interactive HTML dotplots from MAG-bias attrition lineage CSVs."
    )
    ap.add_argument("--results_dir", required=True, help="Pipeline output directory (e.g. results_mag_bias/)")
    ap.add_argument("--outname_all", default="interactive_attrition_ALL.html", help="Filename for combined plot")
    ap.add_argument("--min_reads", type=int, default=5, help="Minimum Total Count to include a lineage")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        raise SystemExit(f"[ERROR] results_dir not found: {results_dir}")

    csvs = find_attrition_csvs(results_dir)
    if not csvs:
        raise SystemExit(f"[ERROR] No *_attrition_lineage*.csv found under {results_dir}")

    all_frames = []

    for csv_path in csvs:
        sample_id = csv_path.parent.name
        df = build_dataframe(csv_path, sample_id)
        df = df[df["Total Count"] >= args.min_reads].copy()
        if df.empty:
            continue

        fig = plot_for_one_sample(df, title=f"Attrition dotplot — {sample_id}")
        out_html = csv_path.parent / f"interactive_attrition_{sample_id}.html"
        write_clickable_html_with_slider(fig, out_html)
        print(f"[OK] Wrote {out_html}")

        all_frames.append(df)

    if not all_frames:
        raise SystemExit("[ERROR] All CSVs were empty after filtering. Try lowering --min_reads.")

    df_all = pd.concat(all_frames, ignore_index=True)
    fig_all = plot_for_one_sample(df_all, title="Attrition dotplot — ALL samples")
    out_all = results_dir / args.outname_all
    write_clickable_html_with_slider(fig_all, out_all)
    print(f"[OK] Wrote {out_all}")


if __name__ == "__main__":
    main()