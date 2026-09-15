// Mirror of backend/app/palette.py. Keep in sync.
// dataviz skill validated colorblind-safe categorical hues; grey = annotated.

export const CATEGORY_COLORS: Record<string, { light: string; dark: string }> = {
  annotated: { light: "#898781", dark: "#898781" },
  exon_skipping: { light: "#2a78d6", dark: "#3987e5" },
  alt_5p: { light: "#eb6834", dark: "#d95926" },
  alt_3p: { light: "#1baf7a", dark: "#199e70" },
  isoform_switch: { light: "#eda100", dark: "#c98500" },
  novel_exon: { light: "#008300", dark: "#008300" },
  novel: { light: "#4a3aa7", dark: "#9085e9" },
};

export const CATEGORY_LABELS: Record<string, string> = {
  annotated: "Annotated",
  exon_skipping: "Exon skipping",
  alt_5p: "Alt 5' splice site",
  alt_3p: "Alt 3' splice site",
  isoform_switch: "Isoform switch",
  novel_exon: "Novel exon",
  novel: "Novel / unknown",
};

const prefersDark = () =>
  typeof window !== "undefined" &&
  window.matchMedia &&
  window.matchMedia("(prefers-color-scheme: dark)").matches;

export function categoryColor(category: string): string {
  const c = CATEGORY_COLORS[category] ?? CATEGORY_COLORS.novel;
  return prefersDark() ? c.dark : c.light;
}
