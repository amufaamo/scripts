library(rtracklayer)
library(tidyverse)

# gffファイルから遺伝子アノテーションを抽出する関数
extract_gene_annotations <- function(gff_file, output_file = "gene_annotations.tsv") {
  # gffファイルの読み込み
  gff <- rtracklayer::readGFF(gff_file)

  # アトリビュート列から情報を抽出する関数
  extract_attribute <- function(attributes, attribute_name) {
    sapply(attributes, function(attr_str) {
      if (is.na(attr_str)) {
        return(NA_character_)
      }
      match <- str_match(attr_str, paste0(attribute_name, "=([^;]+)"))
      if (!is.na(match[1,1])) {
        return(match[1, 2])
      } else {
        return(NA_character_)
      }
    })
  }

  # geneフィーチャーの情報を抽出
  gene_features <- gff %>%
    filter(type == "gene") %>%
    mutate(
      gene_id = extract_attribute(attributes, "ID"),
      locus_tag = extract_attribute(attributes, "locus_tag"),
      name = extract_attribute(attributes, "Name"),
      gene_identifier = ifelse(!is.na(locus_tag), locus_tag, ifelse(!is.na(name), name, sub("gene-", "", gene_id)))
    ) %>%
    select(gene_identifier, seqid, start, end, strand)

  # mRNAまたはCDSフィーチャーの情報を抽出
  mrna_cds_features <- gff %>%
    filter(type %in% c("mRNA", "CDS")) %>%
    mutate(
      parent = extract_attribute(attributes, "Parent"),
      product = extract_attribute(attributes, "product"),
      note = extract_attribute(attributes, "Note")
    ) %>%
    select(parent, product, note)

  # gene_idにproduct, noteを紐付ける
  annotation_dict <- mrna_cds_features %>%
    mutate(gene_identifier = sub("rna-|cds-", "", parent)) %>%
    group_by(gene_identifier) %>%
    summarise(
      product = paste(na.omit(unique(product)), collapse = "; "),
      note = paste(na.omit(unique(note)), collapse = "; ")
    )

  # 遺伝子情報にproduct, noteを結合
  gene_annotations <- gene_features %>%
    left_join(annotation_dict, by = c("gene_identifier" = "gene_identifier"))

  # タブ区切りファイルの作成
  write_tsv(gene_annotations %>% select(gene_identifier, product, note), file = output_file)

  print(paste("アノテーション情報を", output_file, "に保存しました。"))
  invisible(gene_annotations %>% select(gene_identifier, product, note)) # 結果を返り値として返す (オプション)
}
