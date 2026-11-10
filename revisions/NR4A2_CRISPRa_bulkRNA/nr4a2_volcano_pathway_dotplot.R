# NR4A2 enhancer CRISPRa (NR4A2_1, NR4A2_2) vs SCRAMBLE: volcano with the genes of the top up-regulated pathway highlighted, and GSEA dot plot of the top enriched and depleted pathways.
# pathway = the top up-regulated MSigDB Hallmark set of the saved GSEA (by NES), or PATHWAY=<Term> (as spelled in the GSEA table); NR4A2 itself is in its leading edge.
# volcano: all genes grey; pathway members (msigdbr Hallmark) pale red; leading-edge genes (GSEA Lead_genes) red and labelled (top N_LAB by significance); NR4A2 called out with its status in the pathway
# dot plot: top N_UP enriched and top N_DOWN depleted Hallmark pathways (FDR < 0.05) by NES; colour = NES (red enriched, blue depleted), size = -log10 FDR (FDR 0 capped at half the smallest non-zero FDR);
#           pathways with NR4A2 in their leading edge are ringed in black and labelled (sets that only contain NR4A2 are not marked)
# inputs : <DE_DIR>/DESeq2_results_full.csv, <DE_DIR>/gsea/gsea_MSigDB_Hallmark_2020.csv
# outputs: <DE_DIR>/05_volcano_top_pathway_and_dotplot.{png,pdf}, 05_top_pathway_genes.csv
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(patchwork) })
DE_DIR <- Sys.getenv("DE_DIR", unset = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/de/NR4A2_1_2_vs_SCRAMBLE")
N_UP <- as.integer(Sys.getenv("N_UP", unset = "10")); N_DOWN <- as.integer(Sys.getenv("N_DOWN", unset = "7")); N_LAB <- as.integer(Sys.getenv("N_LAB", unset = "18")); RED <- "#c8322b"; BLUE <- "#2f6db5"
de <- read.csv(file.path(DE_DIR, "DESeq2_results_full.csv"), stringsAsFactors = FALSE); de <- de[!is.na(de$padj) & !is.na(de$symbol) & de$symbol != "", ]; de$mlp <- -log10(pmax(de$padj, 1e-300))
gs <- read.csv(file.path(DE_DIR, "gsea", "gsea_MSigDB_Hallmark_2020.csv"), check.names = FALSE, stringsAsFactors = FALSE); gs$n_lead <- as.integer(sub("/.*", "", gs$`Tag %`)); gs$FDR <- gs$`FDR q-val`
sig <- gs[gs$FDR < 0.05, ]; PW <- Sys.getenv("PATHWAY", unset = sig$Term[which.max(sig$NES)]); stopifnot(PW %in% gs$Term); lead <- strsplit(gs$Lead_genes[gs$Term == PW], ";")[[1]]
# Hallmark membership from msigdbr (set name matched to the Enrichr term); without it only the leading edge is known
key <- function(x) gsub("[^A-Z0-9]", "", toupper(gsub("alpha", "A", gsub("beta", "B", x, ignore.case = TRUE), ignore.case = TRUE)))
H <- tryCatch({ h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, key(sub("^HALLMARK_", "", h$gs_name))) }, error = function(e) { message("msigdbr not used: ", conditionMessage(e)); list() })
memb <- if (!is.null(H[[key(PW)]])) unique(H[[key(PW)]]) else lead
gs$NR4A2 <- ifelse(sapply(strsplit(gs$Lead_genes, ";"), function(g) "NR4A2" %in% g), "leading edge", ifelse(sapply(key(gs$Term), function(k) "NR4A2" %in% H[[k]]), "member", "no")); NR_IN <- gs$NR4A2[gs$Term == PW]
cat("Hallmark sets containing NR4A2:\n"); print(gs[gs$NR4A2 != "no", c("Term", "NES", "FDR", "NR4A2")], row.names = FALSE)
de$cls <- factor(ifelse(de$symbol %in% lead, "leading edge", ifelse(de$symbol %in% memb, "other pathway genes", "other")), levels = c("other", "other pathway genes", "leading edge"))
pg <- de[de$cls != "other", c("symbol", "cls", "baseMean", "log2FoldChange", "pvalue", "padj")]; pg <- pg[order(pg$cls != "leading edge", pg$padj), ]; write.csv(pg, file.path(DE_DIR, "05_top_pathway_genes.csv"), row.names = FALSE)
le <- de[de$cls == "leading edge", ]; lab <- le[order(le$padj), ][seq_len(min(N_LAB, nrow(le))), ]; lab <- lab[lab$symbol != "NR4A2", ]
cat(sprintf("pathway: %s | NES %.2f, FDR %.3g | %d leading-edge genes (%d in DE table, %d with padj < 0.05 and log2FC > 0) | %d set members in DE table\n", PW, gs$NES[gs$Term == PW], gs$FDR[gs$Term == PW], length(lead), nrow(le), sum(le$padj < 0.05 & le$log2FoldChange > 0), sum(de$cls != "other")))
YCAP <- max(de$mlp[de$cls != "other"]) * 1.12; n_cap <- sum(de$mlp > YCAP); de$mlp <- pmin(de$mlp, YCAP); lab$mlp <- pmin(lab$mlp, YCAP)        # y axis capped just above the most significant pathway gene
xl <- quantile(abs(de$log2FoldChange), 0.9995) * 1.05; de$x <- pmax(pmin(de$log2FoldChange, xl), -xl); lab$x <- pmax(pmin(lab$log2FoldChange, xl), -xl); nr <- de[de$symbol == "NR4A2", ]
th <- theme_classic(base_size = 8.5) + theme(plot.title = element_text(face = "bold", size = 9.5), plot.subtitle = element_text(size = 7, colour = "grey30"), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7), legend.key.size = unit(3.2, "mm"), plot.margin = margin(4, 6, 4, 4))
p1 <- ggplot(de, aes(x, mlp)) + geom_point(data = de[de$cls == "other", ], colour = "grey82", size = 0.35, alpha = 0.6) + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(data = de[de$cls == "other pathway genes", ], aes(fill = cls), shape = 21, colour = "grey40", size = 1.3, stroke = 0.25) + geom_point(data = de[de$cls == "leading edge", ], aes(fill = cls), shape = 21, colour = "grey15", size = 1.9, stroke = 0.25) +
  geom_point(data = nr, shape = 23, fill = "black", colour = "black", size = 2.8) +
  geom_text_repel(data = lab, aes(label = symbol), colour = "grey10", size = 2.4, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40", box.padding = 0.25, seed = 1) +
  geom_label_repel(data = nr, aes(label = sprintf("NR4A2: %s\nlog2FC %+.2f, adj. p %.3f", c(`leading edge` = "leading-edge gene", member = "pathway member", no = "not in pathway")[[NR_IN]], log2FoldChange, padj)), colour = "black", fontface = "bold", size = 2.5, fill = "white", label.size = 0.3,
                   nudge_x = -2.3, nudge_y = YCAP * 0.4, segment.colour = "black", segment.size = 0.4, arrow = arrow(length = unit(1.3, "mm"), type = "closed"), min.segment.length = 0, lineheight = 0.95, seed = 1) +
  scale_fill_manual(values = c(`other pathway genes` = "#f3c4bf", `leading edge` = RED), name = NULL, labels = c(sprintf("other set members (%d)", sum(de$cls == "other pathway genes")), sprintf("leading edge (%d)", nrow(le)))) +
  labs(x = "log2 fold change, NR4A2 enhancer CRISPRa vs scramble", y = "-log10 adjusted p", title = PW, subtitle = sprintf("top enriched pathway (NES %.2f, FDR %s); dashed = adj. p 0.05", gs$NES[gs$Term == PW], ifelse(gs$FDR[gs$Term == PW] == 0, "< 0.001", signif(gs$FDR[gs$Term == PW], 2)))) +
  th + theme(legend.position = c(0.01, 1), legend.justification = c(0, 1), legend.background = element_blank())
sig <- gs[gs$FDR < 0.05, ]; sig$Term <- sub("^Pperoxisome$", "Peroxisome", sig$Term)   # after the NR4A2 column was added; "Pperoxisome" is a typo of the Enrichr library
up <- sig[sig$NES > 0, ]; up <- head(up[order(-up$NES), ], N_UP); dn <- sig[sig$NES < 0, ]; dn <- head(dn[order(dn$NES), ], N_DOWN); top <- rbind(up, dn); top <- top[order(top$NES), ]; top$Term <- factor(top$Term, levels = top$Term)
top$mlf <- -log10(pmax(top$FDR, min(gs$FDR[gs$FDR > 0]) / 2)); nl <- max(abs(top$NES)); tn <- top[top$NR4A2 == "leading edge", ]   # annotate only where NR4A2 is in the leading edge, not sets that merely contain it
p2 <- ggplot(top, aes(NES, Term)) + geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) + geom_segment(aes(x = 0, xend = NES, yend = Term), colour = "grey85", linewidth = 0.3) + geom_point(aes(size = mlf, fill = NES), shape = 21, colour = "grey25", stroke = 0.25) +
  geom_point(data = tn, aes(size = mlf), shape = 21, fill = NA, colour = "black", stroke = 1.1, show.legend = FALSE) + geom_text(data = tn, aes(x = NES + ifelse(NES > 0, 0.22, -0.22), hjust = ifelse(NES > 0, 0, 1), label = "NR4A2"), size = 2.4, fontface = "bold") +
  scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0, limits = c(-nl, nl), name = "NES") + scale_size_continuous(range = c(1.6, 5), name = "-log10 FDR") + scale_x_continuous(limits = c(-nl - 0.5, nl + 1.3)) + coord_cartesian(clip = "off") +
  labs(x = "normalised enrichment score", y = NULL, title = "GSEA, Hallmark", subtitle = "red = enriched, blue = depleted; ring = NR4A2 in leading edge") + th + theme(axis.text.y = element_text(size = 7.5, colour = "grey10"), legend.box = "vertical", legend.spacing.y = unit(1, "mm"))
fig <- p1 + p2 + plot_layout(widths = c(1.15, 1)); for (ext in c("png", "pdf")) ggsave(file.path(DE_DIR, paste0("05_volcano_top_pathway_and_dotplot.", ext)), fig, width = 9, height = 4.2, dpi = 300)
cat(sprintf("saved 05_volcano_top_pathway_and_dotplot (%d enriched, %d depleted pathways shown; %d volcano genes above the y limit)\n", nrow(up), nrow(dn), n_cap))
