library(data.table)
library(stringr)
library(R.utils) 

input_folder <- "/Users/awsms1/RNA lab/aligned clusters/400_aln"  
output_folder <- "/Users/awsms1/RNA lab/Cleaned aligned clusters/400"  

if (!dir.exists(output_folder)) dir.create(output_folder)

aln_files <- list.files(input_folder, pattern = "\\.aln$", full.names = TRUE)

# ---- Process each file ----
for (input_file in aln_files) {
  
  output_file <- file.path(output_folder, 
                           paste0(tools::file_path_sans_ext(basename(input_file)), "_cleaned.csv"))
  
  total_lines <- countLines(input_file)
  
  results <- list()
  current_query <- NA
  
  con <- file(input_file, "r")
  pb <- txtProgressBar(min = 0, max = total_lines, style = 3)
  
  line_num <- 0
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break  # EOF
    line_num <- line_num + 1
    
    line <- str_trim(line)
    
    if (startsWith(line, "Query")) {
      current_query <- str_trim(sub("\\s*\\[.*", "", str_replace(line, "^Query\\s+>", "")))
      
    } else if (grepl("^\\d+%\\s+", line)) {
      
      parts <- str_split_fixed(line, "\\s+", 3)
      
      if (nchar(parts[3]) > 0) {
        target_clean <- str_trim(sub("\\s*\\[.*", "", parts[3]))
        
        results[[length(results) + 1]] <- list(
          Query_ID = current_query,
          PercentId = parts[1],
          TLen = parts[2],
          Target = target_clean
        )
      }
    }
    
    if (line_num %% 1000 == 0 || line_num == total_lines) {
      setTxtProgressBar(pb, line_num)
    }
  }
  
  close(con)
  close(pb)
  
  df <- rbindlist(results)
  
  fwrite(df, output_file)
  cat("\nProcessed file saved to:", output_file, "\n")
}
