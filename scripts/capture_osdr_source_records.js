#!/usr/bin/env node
/**
 * Capture official OSDR source-record pages for BioVis proof panels.
 *
 * The output is intentionally treated as source-record evidence only. It does
 * not validate local matrices, splits, or benchmark performance.
 */

const fs = require("fs");
const path = require("path");
const { chromium } = require("playwright");

const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(ROOT, "output", "biovis_osdr_source_record_captures_v0_1");
const CREATED = "2026-06-02";

const SOURCES = [
  {
    source_id: "OSD-863",
    lane: "human_organoid",
    url: "https://osdr.nasa.gov/bio/repo/data/studies/OSD-863",
    expected_terms: ["OSD-863", "GLDS-716", "GSE259421", "organoid"],
  },
  {
    source_id: "OSD-871",
    lane: "human_organoid",
    url: "https://osdr.nasa.gov/bio/repo/data/studies/OSD-871",
    expected_terms: ["OSD-871", "GLDS-720", "GSE259421", "organoid"],
  },
  {
    source_id: "OSD-120",
    lane: "multispecies",
    url: "https://osdr.nasa.gov/bio/repo/data/studies/OSD-120",
    expected_terms: ["OSD-120", "GLDS-120", "Arabidopsis"],
  },
  {
    source_id: "OSD-918",
    lane: "single_cell",
    url: "https://osdr.nasa.gov/bio/repo/data/studies/OSD-918",
    expected_terms: ["OSD-918", "RRRM", "single"],
  },
];

function rel(filePath) {
  return path.relative(ROOT, filePath);
}

function ensureDir(dir) {
  fs.mkdirSync(dir, { recursive: true });
}

function safeName(sourceId, suffix) {
  return `${sourceId.toLowerCase().replace(/[^a-z0-9]+/g, "_")}_${suffix}`;
}

function localBrowserExecutable() {
  const candidates = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
  ];
  return candidates.find((candidate) => fs.existsSync(candidate)) || "";
}

async function captureSource(browser, source) {
  const page = await browser.newPage({
    viewport: { width: 1600, height: 1200 },
    deviceScaleFactor: 1,
  });
  const record = {
    source_id: source.source_id,
    lane: source.lane,
    url: source.url,
    expected_terms: source.expected_terms,
    captured_at: new Date().toISOString(),
    status: "pending",
    title: "",
    visible_text_sample: "",
    expected_term_hits: {},
    viewport_screenshot: "",
    fullpage_screenshot: "",
    notes: "",
  };
  try {
    await page.goto(source.url, { waitUntil: "domcontentloaded", timeout: 45000 });
    await page.waitForTimeout(5000);
    const title = await page.title();
    const text = await page.locator("body").innerText({ timeout: 15000 }).catch(() => "");
    record.title = title;
    record.visible_text_sample = text.slice(0, 4000);
    for (const term of source.expected_terms) {
      record.expected_term_hits[term] = text.toLowerCase().includes(term.toLowerCase()) || title.toLowerCase().includes(term.toLowerCase());
    }
    const viewport = path.join(OUT, safeName(source.source_id, "viewport.png"));
    const fullpage = path.join(OUT, safeName(source.source_id, "fullpage.png"));
    await page.screenshot({ path: viewport, fullPage: false });
    await page.screenshot({ path: fullpage, fullPage: true });
    record.viewport_screenshot = rel(viewport);
    record.fullpage_screenshot = rel(fullpage);
    record.status = Object.values(record.expected_term_hits).some(Boolean) ? "captured_review_required" : "captured_term_miss_review_required";
  } catch (error) {
    record.status = "capture_failed";
    record.notes = String(error && error.message ? error.message : error);
  } finally {
    await page.close();
  }
  return record;
}

async function main() {
  ensureDir(OUT);
  const executablePath = localBrowserExecutable();
  const browser = await chromium.launch({
    headless: true,
    ...(executablePath ? { executablePath } : {}),
  });
  const records = [];
  try {
    for (const source of SOURCES) {
      records.push(await captureSource(browser, source));
    }
  } finally {
    await browser.close();
  }
  const manifest = {
    created: CREATED,
    scope: "Official OSDR source-record screenshots for v0.4 proof replacement.",
    records,
    claim_boundary:
      "Screenshots prove source-record page availability and accession context only; they do not validate local payload completeness or benchmark claims.",
  };
  fs.writeFileSync(path.join(OUT, "osdr_source_record_capture_manifest.json"), `${JSON.stringify(manifest, null, 2)}\n`);
  const summary = records.map((record) => ({
    source_id: record.source_id,
    status: record.status,
    viewport_screenshot: record.viewport_screenshot,
    fullpage_screenshot: record.fullpage_screenshot,
    hit_count: Object.values(record.expected_term_hits).filter(Boolean).length,
  }));
  fs.writeFileSync(path.join(OUT, "osdr_source_record_capture_summary.json"), `${JSON.stringify(summary, null, 2)}\n`);
  console.log(JSON.stringify({ output: rel(OUT), records: records.length, statuses: summary }, null, 2));
}

main().catch((error) => {
  console.error(error);
  process.exit(1);
});
