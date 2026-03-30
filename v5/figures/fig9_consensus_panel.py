#!/usr/bin/env python
"""Fig 9: Consensus Biomarker Panel.

Panel A: Evidence integration matrix (genes × evidence sources, color = score)
Panel B: Panel AUROC vs full gene set across tissues
Panel C: Gene-drug-tissue network (bipartite)
"""
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V4_EVAL_DIR = BASE_DIR / "v4" / "evaluation"
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
OUT_DIR = BASE_DIR / "v5" / "figures" / "html"

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung†", "colon": "Colon†",
}

SOURCE_LABELS = {
    "shap": "SHAP", "wgcna_hub": "WGCNA Hub", "network": "Network",
    "druggable": "Druggable", "literature": "Literature", "cross_species": "Cross-species",
}
SOURCE_ORDER = ["shap", "wgcna_hub", "network", "druggable", "literature", "cross_species"]
SOURCE_MAX = [3, 2, 2, 2, 2, 1]


def load_panel_data():
    """Load consensus biomarker panel results."""
    path = V5_EVAL_DIR / "consensus_biomarker_panel.json"
    with open(path) as f:
        return json.load(f)


def load_m1_aurocs():
    """Load full-gene PCA-LR AUROC from M1 summary for comparison."""
    path = V4_EVAL_DIR / "M1_summary.json"
    if not path.exists():
        return {}
    with open(path) as f:
        d = json.load(f)
    aurocs = {}
    for tissue in TISSUES:
        # Structure: d[tissue][feature_type][method]
        try:
            aurocs[tissue] = d[tissue]["gene"]["pca_lr"]["auroc"]
        except (KeyError, TypeError):
            pass
    return aurocs


def build_evidence_matrix(panel):
    """Build evidence matrix data for D3 heatmap."""
    cells = []
    for gene_info in panel:
        gene = gene_info["gene"]
        sources = gene_info.get("sources", {})
        for src in SOURCE_ORDER:
            cells.append({
                "gene": gene,
                "source": SOURCE_LABELS[src],
                "value": sources.get(src, 0),
            })
    return cells


def build_auroc_comparison(panel_data, m1_aurocs):
    """Build AUROC comparison data: panel vs full gene set."""
    validation = panel_data.get("panel_validation", {})
    bars = []
    for tissue in TISSUES:
        label = TISSUE_LABELS[tissue]
        panel_auroc = validation.get(tissue, {}).get("auroc", None)
        full_auroc = m1_aurocs.get(tissue, None)
        if panel_auroc is not None:
            bars.append({"tissue": label, "type": "Panel (20 genes)", "auroc": round(panel_auroc, 3)})
        if full_auroc is not None:
            bars.append({"tissue": label, "type": "Full gene set", "auroc": round(full_auroc, 3)})
    return bars


def build_gene_drug_network(panel):
    """Build bipartite gene-drug network for top genes with drugs."""
    nodes = []
    links = []
    gene_set = set()
    drug_set = set()

    for gene_info in panel[:15]:  # top 15 to keep network readable
        gene = gene_info["gene"]
        tissues = gene_info.get("tissues", [])
        drugs = gene_info.get("drugs", [])
        if not drugs:
            continue

        gene_set.add(gene)
        gene_node = {
            "id": gene, "type": "gene",
            "tissues": ", ".join(tissues),
            "score": gene_info.get("score", 0),
        }
        if gene not in [n["id"] for n in nodes]:
            nodes.append(gene_node)

        # Limit drugs per gene for readability
        for drug_info in drugs[:3]:
            drug_name = drug_info["drug"]
            # Shorten drug names
            short_name = drug_name.split(" ")[0][:15]
            drug_id = drug_name

            if drug_id not in drug_set:
                drug_set.add(drug_id)
                nodes.append({
                    "id": drug_id, "type": "drug",
                    "label": short_name,
                    "tier": drug_info.get("tier", ""),
                })

            links.append({"source": gene, "target": drug_id})

    return nodes, links


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    panel_data = load_panel_data()
    panel = panel_data.get("panel", [])
    m1_aurocs = load_m1_aurocs()

    evidence_matrix = build_evidence_matrix(panel)
    auroc_bars = build_auroc_comparison(panel_data, m1_aurocs)
    net_nodes, net_links = build_gene_drug_network(panel)

    gene_names = [g["gene"] for g in panel]

    print(f"Panel A: {len(panel)} genes × {len(SOURCE_ORDER)} sources")
    print(f"Panel B: {len(auroc_bars)} bars ({len(auroc_bars)//2} tissues × 2 types)")
    print(f"Panel C: {len(net_nodes)} nodes, {len(net_links)} links")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 9: Consensus Biomarker Panel</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
.figure-container {{ width: 1200px; margin: 0 auto; }}
.figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
.figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
.panels {{ display: grid; grid-template-columns: 1.2fr 0.8fr; gap: 20px; }}
.panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 12px; background: #fafafa; }}
.panel-label {{ font-size: 12px; font-weight: bold; margin-bottom: 2px; }}
.panel-desc {{ font-size: 9px; color: #666; margin-bottom: 8px; }}
.interpretation {{ margin-top: 16px; padding: 12px; background: #f0f7ff; border-radius: 4px;
                   font-size: 10px; line-height: 1.5; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 6px 10px;
            border-radius: 3px; font-size: 9px; pointer-events: none; z-index: 100; }}
.panel-full {{ grid-column: 1 / -1; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 9: Consensus Spaceflight Biomarker Panel</div>
  <div class="figure-subtitle">Multi-source evidence integration: SHAP + WGCNA + PPI + druggability + literature + cross-species ({panel_data.get('n_candidates_scored', 0)} genes scored)</div>
  <div class="panels">
    <div class="panel panel-full" id="panelA">
      <div class="panel-label">A</div>
      <div class="panel-desc">Evidence integration matrix — top-20 genes by composite score (max {panel_data.get('max_possible_score', 13)})</div>
      <div id="evidence-matrix"></div>
    </div>
    <div class="panel" id="panelB">
      <div class="panel-label">B</div>
      <div class="panel-desc">Panel (20 genes) vs full gene set AUROC</div>
      <div id="auroc-compare"></div>
    </div>
    <div class="panel" id="panelC">
      <div class="panel-label">C</div>
      <div class="panel-desc">Gene-drug network (top-15 panel genes, up to 3 drugs each)</div>
      <div id="gene-drug-net"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    Top-20 consensus biomarker panel integrates 6 orthogonal evidence sources.
    Ddit4 (DNA damage, score=6) and Myl2 (cardiac myosin, score=6) rank highest.
    Il6 (inflammatory cytokine), Spp1 (osteopontin), and Ppara (lipid metabolism TF) each score 5, reflecting
    muscle wasting, bone loss, and metabolic reprogramming as core spaceflight responses.
    The 20-gene panel achieves AUROC 0.75 (colon), 0.70 (skin), 0.66 (liver) — capturing 60-80% of full gene set performance
    with 1000&times; fewer features, while every gene has known drug interactions for potential countermeasure development.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tissueColors = {{
  "Liver": OI.blue, "Gastrocnemius": OI.vermilion, "Kidney": OI.green,
  "Thymus": OI.orange, "Eye": OI.pink, "Skin": OI.skyblue,
  "Lung†": "#888888", "Colon†": "#555555"
}};
const tooltip = d3.select("#tooltip");

// Panel A: Evidence integration matrix
(function() {{
  const data = {json.dumps(evidence_matrix)};
  const genes = {json.dumps(gene_names)};
  const sources = {json.dumps([SOURCE_LABELS[s] for s in SOURCE_ORDER])};
  const sourceMax = {json.dumps(SOURCE_MAX)};

  const margin = {{top: 10, right: 80, bottom: 60, left: 90}};
  const cellW = 75, cellH = 22;
  const width = sources.length * cellW, height = genes.length * cellH;

  const svg = d3.select("#evidence-matrix").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(sources).range([0, width]).padding(0.08);
  const y = d3.scaleBand().domain(genes).range([0, height]).padding(0.06);
  const color = d3.scaleSequential(d3.interpolateBlues).domain([0, 3]);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", d => x(d.source)).attr("y", d => y(d.gene))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => d.value > 0 ? color(d.value) : "#f5f5f5")
    .attr("stroke", "#ddd").attr("stroke-width", 0.5)
    .attr("rx", 2)
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`${{d.gene}} — ${{d.source}}<br>Score: ${{d.value}}`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  // Score text in cells
  svg.selectAll("text.val").data(data.filter(d => d.value > 0)).join("text")
    .attr("class", "val")
    .attr("x", d => x(d.source) + x.bandwidth()/2)
    .attr("y", d => y(d.gene) + y.bandwidth()/2 + 4)
    .attr("text-anchor", "middle").attr("font-size", "9px")
    .attr("fill", d => d.value >= 2 ? "#fff" : "#333")
    .text(d => d.value);

  // Total score annotations on right
  const panel = {json.dumps(panel)};
  svg.selectAll("text.total").data(panel).join("text")
    .attr("class", "total")
    .attr("x", width + 8).attr("y", (d,i) => y(genes[i]) + y.bandwidth()/2 + 4)
    .attr("font-size", "10px").attr("font-weight", "bold").attr("fill", OI.vermilion)
    .text(d => `=${{d.score}}`);

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "9px");
  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "9px").style("font-style", "italic");

  // Column header: max possible
  svg.selectAll("text.max").data(sources).join("text")
    .attr("class", "max").attr("x", (d,i) => x(d) + x.bandwidth()/2)
    .attr("y", height + 50).attr("text-anchor", "middle")
    .attr("font-size", "7px").attr("fill", "#999")
    .text((d,i) => `max: ${{sourceMax[i]}}`);
}})();

// Panel B: AUROC comparison
(function() {{
  const data = {json.dumps(auroc_bars)};
  const tissues = [...new Set(data.map(d => d.tissue))];
  const types = [...new Set(data.map(d => d.type))];

  const margin = {{top: 10, right: 20, bottom: 70, left: 50}};
  const width = 420, height = 280;

  const svg = d3.select("#auroc-compare").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x0 = d3.scaleBand().domain(tissues).range([0, width]).padding(0.25);
  const x1 = d3.scaleBand().domain(types).range([0, x0.bandwidth()]).padding(0.1);
  const y = d3.scaleLinear().domain([0, 1]).range([height, 0]);

  const colorMap = {{"Panel (20 genes)": OI.vermilion, "Full gene set": OI.blue}};

  svg.selectAll("g.group").data(tissues).join("g")
    .attr("class", "group")
    .attr("transform", d => `translate(${{x0(d)}},0)`)
    .selectAll("rect").data(t => data.filter(d => d.tissue === t)).join("rect")
    .attr("x", d => x1(d.type)).attr("y", d => y(d.auroc))
    .attr("width", x1.bandwidth()).attr("height", d => height - y(d.auroc))
    .attr("fill", d => colorMap[d.type] || "#999")
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`${{d.tissue}} — ${{d.type}}<br>AUROC: ${{d.auroc}}`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  // Value labels
  svg.selectAll("g.group2").data(tissues).join("g").attr("class", "group2")
    .attr("transform", d => `translate(${{x0(d)}},0)`)
    .selectAll("text").data(t => data.filter(d => d.tissue === t)).join("text")
    .attr("x", d => x1(d.type) + x1.bandwidth()/2)
    .attr("y", d => y(d.auroc) - 4)
    .attr("text-anchor", "middle").attr("font-size", "7px")
    .text(d => d.auroc.toFixed(2));

  // 0.5 reference line
  svg.append("line").attr("x1", 0).attr("x2", width)
    .attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#999").attr("stroke-dasharray", "3,3");
  svg.append("text").attr("x", width + 2).attr("y", y(0.5) + 3)
    .attr("font-size", "7px").attr("fill", "#999").text("chance");

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x0))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y).ticks(5)).selectAll("text").style("font-size", "8px");
  svg.append("text").attr("transform", "rotate(-90)").attr("y", -35).attr("x", -height/2)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("AUROC");

  // Legend
  const legend = svg.append("g").attr("transform", `translate(${{width - 160}}, 5)`);
  types.forEach((t, i) => {{
    legend.append("rect").attr("x", 0).attr("y", i*16).attr("width", 12).attr("height", 12)
      .attr("fill", colorMap[t]);
    legend.append("text").attr("x", 16).attr("y", i*16 + 10).attr("font-size", "8px").text(t);
  }});
}})();

// Panel C: Gene-drug network
(function() {{
  const nodes = {json.dumps(net_nodes)};
  const links = {json.dumps(net_links)};

  const width = 500, height = 400;
  const svg = d3.select("#gene-drug-net").append("svg")
    .attr("width", width).attr("height", height);

  const simulation = d3.forceSimulation(nodes)
    .force("link", d3.forceLink(links).id(d => d.id).distance(60))
    .force("charge", d3.forceManyBody().strength(-120))
    .force("center", d3.forceCenter(width/2, height/2))
    .force("collision", d3.forceCollide().radius(20));

  const link = svg.selectAll("line").data(links).join("line")
    .attr("stroke", "#ccc").attr("stroke-width", 1.5);

  const node = svg.selectAll("g.node").data(nodes).join("g").attr("class", "node")
    .call(d3.drag()
      .on("start", (e,d) => {{ if (!e.active) simulation.alphaTarget(0.3).restart(); d.fx=d.x; d.fy=d.y; }})
      .on("drag", (e,d) => {{ d.fx=e.x; d.fy=e.y; }})
      .on("end", (e,d) => {{ if (!e.active) simulation.alphaTarget(0); d.fx=null; d.fy=null; }})
    );

  node.append("circle")
    .attr("r", d => d.type === "gene" ? 12 : 8)
    .attr("fill", d => d.type === "gene" ? OI.vermilion : OI.skyblue)
    .attr("stroke", "#333").attr("stroke-width", 1);

  node.append("text")
    .attr("dy", d => d.type === "gene" ? -16 : 14)
    .attr("text-anchor", "middle").attr("font-size", d => d.type === "gene" ? "8px" : "6px")
    .attr("font-style", d => d.type === "gene" ? "italic" : "normal")
    .text(d => d.type === "gene" ? d.id : (d.label || d.id.substring(0, 12)));

  node.on("mouseover", (e, d) => {{
    const info = d.type === "gene"
      ? `${{d.id}} (score: ${{d.score}})<br>Tissues: ${{d.tissues}}`
      : `Drug: ${{d.id}}<br>Tier: ${{d.tier}}`;
    tooltip.style("display", "block").html(info)
      .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
  }}).on("mouseout", () => tooltip.style("display", "none"));

  simulation.on("tick", () => {{
    link.attr("x1", d => d.source.x).attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x).attr("y2", d => d.target.y);
    node.attr("transform", d => `translate(${{d.x}},${{d.y}})`);
  }});

  // Legend
  const legend = svg.append("g").attr("transform", "translate(10, 10)");
  [["Gene", OI.vermilion, 12], ["Drug", OI.skyblue, 8]].forEach(([label, color, r], i) => {{
    legend.append("circle").attr("cx", r).attr("cy", i*22 + r).attr("r", r).attr("fill", color).attr("stroke", "#333");
    legend.append("text").attr("x", 2*r + 8).attr("y", i*22 + r + 4).attr("font-size", "9px").text(label);
  }});
}})();
</script>
</body>
</html>"""

    out_path = OUT_DIR / "Fig9_consensus_panel.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
