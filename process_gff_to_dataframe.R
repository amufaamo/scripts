# GFFファイルを解析し、geneとCDS情報を統合してdata.frameを返す単一関数
process_gff_to_dataframe <- function(gff_file) {

  # --- 1. GFFファイルの解析とデータの格納 ---
  gene_data <- list()
  cds_data <- list()
  all_gene_keys <- character()
  all_cds_keys <- character()

  if (!file.exists(gff_file)) {
    stop(paste("Error: Input GFF file not found:", gff_file))
  }

  lines <- readLines(gff_file)
  line_num <- 0
  for (line in lines) {
    line_num <- line_num + 1
    line <- trimws(line)
    if (nchar(line) == 0 || substr(line, 1, 1) == "#") { next }
    parts <- strsplit(line, '\t')[[1]]
    if (length(parts) < 9) {
      warning(paste("Skipping line", line_num, "due to insufficient columns:", line))
      next
    }

    feature_type <- parts[3]
    attributes_str <- parts[9]
    attributes <- list()
    current_keys <- character()

    attr_pairs <- strsplit(attributes_str, ';')[[1]]
    for (attr_pair_str in attr_pairs) {
      attr_pair_str <- trimws(attr_pair_str)
      if (nchar(attr_pair_str) == 0) { next }
      attr_split <- strsplit(attr_pair_str, '=', fixed = TRUE)[[1]]
      if (length(attr_split) >= 2) {
        key <- trimws(attr_split[1])
        value <- trimws(paste(attr_split[-1], collapse = '='))
        if (key %in% names(attributes)) {
          if (!is.list(attributes[[key]])) {
            attributes[[key]] <- list(attributes[[key]])
          }
          attributes[[key]] <- c(attributes[[key]], value)
        } else {
          attributes[[key]] <- value
        }
        current_keys <- c(current_keys, key)
      } else if (length(attr_split) == 1 && nchar(attr_split[1]) > 0) {
          key <- trimws(attr_split[1])
          value <- TRUE # キーのみの場合は TRUE
          attributes[[key]] <- value
          current_keys <- c(current_keys, key)
      }
    }

    if (feature_type == "gene") {
      if ('ID' %in% names(attributes)) {
        gene_id <- attributes[['ID']][1]
        if (gene_id %in% names(gene_data)) {
            warning(paste("Duplicate gene ID found:", gene_id, "on line", line_num, ". Merging attributes, later entries might overwrite."))
            for (key in names(attributes)) { gene_data[[gene_id]][[key]] <- attributes[[key]] }
        } else {
             gene_data[[gene_id]] <- attributes
        }
        all_gene_keys <- union(all_gene_keys, current_keys)
      } else {
          warning(paste("Skipping gene feature on line", line_num, "due to missing ID:", line))
      }
    } else if (feature_type == "CDS") {
      if ('Parent' %in% names(attributes)) {
        parent_id <- attributes[['Parent']][1]
        if (parent_id %in% names(cds_data)) {
            existing_cds_attrs <- cds_data[[parent_id]]
            for (key in names(attributes)) {
                if (key %in% names(existing_cds_attrs)) {
                    if (!is.list(existing_cds_attrs[[key]])) { existing_cds_attrs[[key]] <- list(existing_cds_attrs[[key]]) }
                    new_value <- attributes[[key]]
                    if(is.list(new_value)){ existing_cds_attrs[[key]] <- c(existing_cds_attrs[[key]], unlist(new_value)) }
                    else { existing_cds_attrs[[key]] <- c(existing_cds_attrs[[key]], new_value) }
                } else {
                    existing_cds_attrs[[key]] <- attributes[[key]]
                }
            }
             cds_data[[parent_id]] <- existing_cds_attrs
        } else {
             cds_data[[parent_id]] <- attributes
        }
        all_cds_keys <- union(all_cds_keys, current_keys)
      } else {
          warning(paste("Skipping CDS feature on line", line_num, "due to missing Parent attribute:", line))
      }
    }
  }

  # --- 2. gene情報とCDS情報のマージ ---
  merged_gene_info <- list()
  final_all_keys <- all_gene_keys

  cds_prefixed_keys <- character()
  if(length(all_cds_keys) > 0) {
      keys_to_prefix <- setdiff(all_cds_keys, c("ID", "Parent"))
      cds_prefixed_keys <- paste0("CDS_", keys_to_prefix)
      if ("ID" %in% all_cds_keys) cds_prefixed_keys <- c(cds_prefixed_keys, "CDS_ID")
  }
  final_all_keys <- union(final_all_keys, cds_prefixed_keys)

  for (gene_id in names(gene_data)) {
    merged_info <- gene_data[[gene_id]]
    if (gene_id %in% names(cds_data)) {
      current_cds_info <- cds_data[[gene_id]]
      for (key in names(current_cds_info)) {
        if (key == "Parent") { next }
        if (key == "ID") { new_key <- "CDS_ID" }
        else { new_key <- paste0("CDS_", key) }
        value <- current_cds_info[[key]]
        if (new_key %in% names(merged_info)) {
             if (!is.list(merged_info[[new_key]])) { merged_info[[new_key]] <- list(merged_info[[new_key]]) }
             if (is.list(value)) { merged_info[[new_key]] <- c(merged_info[[new_key]], unlist(value)) }
             else { merged_info[[new_key]] <- c(merged_info[[new_key]], value) }
        } else {
             merged_info[[new_key]] <- value
        }
      }
    }
    merged_gene_info[[gene_id]] <- merged_info
  }

  final_all_keys_sorted <- sort(setdiff(final_all_keys, "ID"))

  # --- 3. data.frameの構築 ---
  if (length(merged_gene_info) == 0) {
    warning("No gene features found or processed. Returning empty data.frame.")
    # 空のデータフレームを定義（ヘッダーだけは設定する）
    df_columns <- c("Gene_ID", final_all_keys_sorted)
    empty_df <- data.frame(matrix(ncol = length(df_columns), nrow = 0))
    colnames(empty_df) <- df_columns
    return(empty_df)
  }

  # 列ごとのベクトルを格納するリストを初期化
  output_df_list <- list()
  output_df_list[["Gene_ID"]] <- names(merged_gene_info)

  # 各属性キーについて、全遺伝子の値を取得してベクトル化
  for (col_key in final_all_keys_sorted) {
    # sapplyを使って各gene_idに対応する値を取得・整形
    column_vector <- sapply(names(merged_gene_info), function(gene_id) {
      info <- merged_gene_info[[gene_id]]
      value_str <- NA # デフォルトはNA

      if (col_key %in% names(info)) {
        value <- info[[col_key]]
        if (is.list(value)) {
          # リスト内の非NULL/非空文字要素を抽出し、";"で結合
          valid_values <- unlist(value)[sapply(unlist(value), function(x) !is.null(x) && nchar(trimws(as.character(x))) > 0)]
          if (length(valid_values) > 0) {
            value_str <- paste(valid_values, collapse = ";")
          } else {
            value_str <- "" # 有効な値がない場合は空文字にする（NAでも良い）
          }
        } else if (!is.null(value) && !is.na(value) && nchar(trimws(as.character(value))) > 0) {
          # 単一の値で、NULL/NA/空文字でない場合
          value_str <- as.character(value) # 文字列として格納
        } else if (is.logical(value)) {
            value_str <- as.character(value) # TRUE/FALSEを文字列に
        } else {
           value_str <- "" # 上記以外（NULLや空文字）の場合は空文字にする
        }
      }
      return(value_str)
    }, USE.NAMES = FALSE) # USE.NAMES=FALSEでベクトルとして返す

    output_df_list[[col_key]] <- column_vector
  }

  # リストからデータフレームを作成 (文字列として扱う)
  output_df <- data.frame(output_df_list, stringsAsFactors = FALSE, check.names = FALSE)

  return(output_df)
}
all_header_keys <- result$all_attribute_keys
output_header <- c("ID", all_header_keys)
write_gene_info_to_tsv_all_attributes(gene_information, ts_output_file, output_header)

cat(paste("Gene information with all attributes written to", ts_output_file, "\n"))
