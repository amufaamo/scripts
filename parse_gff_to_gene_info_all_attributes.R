parse_gff_to_gene_info_all_attributes <- function(gff_file) {
  gene_info <- list()
  all_attribute_keys <- character()

  lines <- readLines(gff_file)
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0 || substr(line, 1, 1) == "#") { # 空行とコメント行をスキップ
      next
    }
    parts <- strsplit(line, '\t')[[1]]
    if (length(parts) < 9) {
      next
    }
    attributes_str <- parts[9]
    attributes <- list()
    attr_pairs <- strsplit(attributes_str, ';')[[1]]
    for (attr_pair_str in attr_pairs) {
      attr_pair_str <- trimws(attr_pair_str)
      if (nchar(attr_pair_str) == 0) {
        next
      }
      attr_split <- strsplit(attr_pair_str, '=')[[1]]
      if (length(attr_split) == 2) {
        key <- attr_split[1]
        value <- attr_split[2]
        
        if (key %in% names(attributes)) {
          if (is.list(attributes[[key]])) {
            attributes[[key]] <- c(attributes[[key]], value)
          } else {
            attributes[[key]] <- list(attributes[[key]], value)
          }
        } else {
          attributes[[key]] <- value
        }
        all_attribute_keys <- union(all_attribute_keys, key)
      }
    }
    
    if ('ID' %in% names(attributes)) {
      gene_id <- attributes[['ID']]
      if (!(gene_id %in% names(gene_info))) {
        gene_info[[gene_id]] <- list()
      }
      gene_info[[gene_id]] <- append(gene_info[[gene_id]], attributes)
    }
  }

  return(list(gene_info = gene_info, all_attribute_keys = sort(all_attribute_keys)))
}

write_gene_info_to_tsv_all_attributes <- function(gene_info, output_file, header) {
  output_conn <- file(output_file, "w")
  writeLines(paste(header, collapse = "\t"), output_conn)
  for (gene_id in names(gene_info)) {
    info <- gene_info[[gene_id]]
    row <- gene_id
    for (col in header[-1]) {
      if (col %in% names(info)) {
        value <- info[[col]]
        if (is.list(value)) {
          row <- c(row, paste(value, collapse = ";"))
        } else {
          row <- c(row, value)
        }
      } else {
        row <- c(row, '')
      }
    }
    writeLines(paste(row, collapse = "\t"), output_conn)
  }
  close(output_conn)
}

# メイン処理
gff_file <- "your_input.gff"
ts_output_file <- "gene_info_all_attributes.tsv"

result <- parse_gff_to_gene_info_all_attributes(gff_file)
gene_information <- result$gene_info
all_header_keys <- result$all_attribute_keys
output_header <- c("ID", all_header_keys)
write_gene_info_to_tsv_all_attributes(gene_information, ts_output_file, output_header)

cat(paste("Gene information with all attributes written to", ts_output_file, "\n"))
