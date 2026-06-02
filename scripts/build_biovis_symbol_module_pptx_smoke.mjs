#!/usr/bin/env node
/**
 * Build a minimal PPTX smoke test for the biological symbol/module SVG assets.
 *
 * This is not a production deck builder. It verifies that representative SVG
 * assets can be embedded into a PPTX and rendered/exported by local tooling.
 */

import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";

const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
const OUT = path.join(ROOT, "output", "biovis_symbol_module_bridge_application", "pptx_smoke");
const ASSET = path.join(ROOT, "assets", "biovis_symbol_module_pack_v0_3");
const NODE_MODULES =
  process.env.CODEX_NODE_MODULES ||
  "/Users/jak4013/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules";

const require = createRequire(import.meta.url);
const pptxgen = require(path.join(NODE_MODULES, "pptxgenjs"));

function readSvgDataUri(filePath) {
  const svg = fs.readFileSync(filePath, "utf8");
  return `data:image/svg+xml;base64,${Buffer.from(svg).toString("base64")}`;
}

function ensure() {
  fs.mkdirSync(OUT, { recursive: true });
}

function addLabel(slide, text, x, y, opts = {}) {
  slide.addText(text, {
    x,
    y,
    w: opts.w || 4.2,
    h: opts.h || 0.26,
    fontFace: "Arial",
    fontSize: opts.size || 9,
    bold: opts.bold || false,
    color: opts.color || "17212B",
    margin: 0,
    breakLine: false,
    fit: "shrink",
  });
}

async function main() {
  ensure();
  const pptx = new pptxgen();
  pptx.layout = "LAYOUT_WIDE";
  pptx.author = "SpaceBio-Bench";
  pptx.subject = "Biological symbol/module SVG editability smoke test";
  pptx.title = "Biological symbol module PPTX smoke test";
  pptx.company = "SpaceBio-Bench";
  pptx.lang = "en-US";
  pptx.theme = {
    headFontFace: "Arial",
    bodyFontFace: "Arial",
    lang: "en-US",
  };

  const slide = pptx.addSlide();
  slide.background = { color: "F7F4EC" };
  slide.addText("SVG editability smoke test", {
    x: 0.45,
    y: 0.28,
    w: 12.4,
    h: 0.45,
    fontFace: "Arial",
    fontSize: 26,
    bold: true,
    color: "17212B",
    margin: 0,
  });
  slide.addText("Representative v0.3 biological symbols and modules embedded as SVG media.", {
    x: 0.45,
    y: 0.78,
    w: 12.0,
    h: 0.28,
    fontFace: "Arial",
    fontSize: 11,
    color: "5D6978",
    margin: 0,
  });

  const moduleSvg = path.join(ASSET, "modules", "svg", "sample_to_feature_stack.svg");
  const darkModuleSvg = path.join(ASSET, "modules_dark", "svg", "space_bio_assay_lane.svg");
  const claimSvg = path.join(ASSET, "modules", "svg", "claim_boundary_bar.svg");
  const iconSvg = path.join(ASSET, "symbols", "svg", "organoid_rosette.svg");
  const badgeSvg = path.join(ASSET, "badges", "svg", "source_proof.svg");

  slide.addShape(pptx.ShapeType.rect, {
    x: 0.45,
    y: 1.24,
    w: 7.4,
    h: 2.15,
    fill: { color: "FBFAF7", transparency: 0 },
    line: { color: "AEB8C5", transparency: 30, width: 1 },
    radius: 0.06,
  });
  slide.addImage({ data: readSvgDataUri(moduleSvg), x: 0.64, y: 1.38, w: 7.02, h: 2.02 });
  addLabel(slide, "Light module SVG", 0.62, 3.54, { bold: true });

  slide.addShape(pptx.ShapeType.rect, {
    x: 8.15,
    y: 1.24,
    w: 4.9,
    h: 2.15,
    fill: { color: "0D1720", transparency: 0 },
    line: { color: "617282", transparency: 10, width: 1 },
    radius: 0.06,
  });
  slide.addImage({ data: readSvgDataUri(darkModuleSvg), x: 8.30, y: 1.43, w: 4.58, h: 1.68 });
  addLabel(slide, "Dark module SVG", 8.22, 3.54, { bold: true });

  slide.addImage({ data: readSvgDataUri(claimSvg), x: 0.58, y: 4.26, w: 5.62, h: 1.24 });
  addLabel(slide, "Claim boundary module", 0.62, 5.62, { bold: true });

  slide.addShape(pptx.ShapeType.rect, {
    x: 7.05,
    y: 4.1,
    w: 2.15,
    h: 1.55,
    fill: { color: "FBFAF7" },
    line: { color: "AEB8C5", transparency: 25, width: 1 },
    radius: 0.06,
  });
  slide.addImage({ data: readSvgDataUri(iconSvg), x: 7.50, y: 4.26, w: 0.88, h: 0.88 });
  addLabel(slide, "Symbol SVG", 7.30, 5.32, { w: 1.6, bold: true });

  slide.addImage({ data: readSvgDataUri(badgeSvg), x: 9.80, y: 4.36, w: 2.05, h: 0.72 });
  addLabel(slide, "Badge SVG", 9.92, 5.32, { w: 1.6, bold: true });

  slide.addText("Smoke criteria: PPTX is created, SVG media are present in ppt/media, and local office tooling can open/export it.", {
    x: 0.45,
    y: 6.72,
    w: 12.4,
    h: 0.32,
    fontFace: "Arial",
    fontSize: 10,
    color: "5D6978",
    margin: 0,
  });

  const outPptx = path.join(OUT, "biovis_symbol_module_svg_editability_smoke.pptx");
  await pptx.writeFile({ fileName: outPptx });
  const manifest = {
    created: "2026-06-02",
    purpose: "PPTX smoke test for SVG embedding of biological symbol/module assets.",
    pptx: path.relative(ROOT, outPptx),
    representative_assets: [
      path.relative(ROOT, moduleSvg),
      path.relative(ROOT, darkModuleSvg),
      path.relative(ROOT, claimSvg),
      path.relative(ROOT, iconSvg),
      path.relative(ROOT, badgeSvg),
    ],
    expected_media_type: "SVG media should remain in ppt/media when embedded by pptxgenjs.",
  };
  fs.writeFileSync(path.join(OUT, "biovis_symbol_module_pptx_smoke_manifest.json"), `${JSON.stringify(manifest, null, 2)}\n`);
  console.log(JSON.stringify({ pptx: path.relative(ROOT, outPptx), assets: manifest.representative_assets.length }, null, 2));
}

await main();
