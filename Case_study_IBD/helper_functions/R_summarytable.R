rst <- function (data, vars = NA, out = NA, file = NA, summ = NA, summ.names = NA, 
          add.median = FALSE, group = NA, group.long = FALSE, group.test = FALSE, 
          group.weights = NA, col.breaks = NA, digits = 2, fixed.digits = FALSE, 
          numformat = formatfunc(digits = digits, big.mark = ""), 
          skip.format = c("notNA(x)", "propNA(x)", "countNA(x)", obs.function), 
          factor.percent = TRUE, factor.counts = TRUE, factor.numeric = FALSE, 
          logical.numeric = FALSE, logical.labels = c("No", "Yes"), 
          labels = NA, title = "Summary Statistics", note = NA, anchor = NA, 
          col.width = NA, col.align = NA, align = NA, note.align = "l", 
          fit.page = "\\textwidth", simple.kable = FALSE, obs.function = NA, 
          stat.cutoffs = c(0.0005, 0.005, 0.05),
          opts = list()) 
{
  list2env(opts, envir = environment())
  if (is.null(colnames(data))) {
    stop("Requires data with variable names or column names.")
  }
  if (!is.na(file) & !is.character(file)) {
    stop("Incorrect file name.")
  }
  if (!is.na(out)) {
    if (!(out %in% c("return", "viewer", "browser", "htmlreturn", 
                     "latex", "kable", "csv", "latexpage"))) {
      stop("Unrecognized option for out.")
    }
  }
  if (!identical(vars, NA) & !is.character(vars)) {
    stop("vars must be a character vector.")
  }
  if (!identical(note, NA) & !is.character(note)) {
    stop("note must be a character vector.")
  }
  if (!identical(anchor, NA) & !is.character(anchor)) {
    stop("anchor must be a character variable.")
  }
  if (min(is.na(col.width)) == 0 & (!is.vector(col.width) | 
                                    !is.numeric(col.width) | sum(is.na(col.width)) > 0)) {
    stop("col.width must be a numeric vector with no missing values.")
  }
  if (min(is.na(col.width)) == 0 & (max(col.width) > 100 | 
                                    min(col.width) < 0)) {
    stop("Elements of col.width must be between 0 and 100.")
  }
  if (!is.logical(add.median)) {
    stop("add.median must be TRUE or FALSE.")
  }
  if (!is.list(summ)) {
    if (min(is.na(summ)) == 0 & (!is.vector(summ) | !is.character(summ) | 
                                 sum(is.na(summ)) > 0)) {
      stop("summ must be a character vector with no missing values.")
    }
  }
  if (!is.list(summ.names)) {
    if (min(is.na(summ.names)) == 0 & (!is.vector(summ.names) | 
                                       !is.character(summ.names) | sum(is.na(summ.names)) > 
                                       0)) {
      stop("summ.names must be a character vector with no missing values.")
    }
  }
  if (!is.na(group) & !is.character(group)) {
    stop("group must be a string referring to a grouping variable in the data.")
    if (!(group %in% colnames(data))) {
      stop("group must be a column name in the data.")
    }
  }
  if (!is.logical(group.test) & !is.list(group.test)) {
    stop("group.test must be TRUE, FALSE, or a named list of options to pass to independence.test's opts argument.")
  }
  if (!identical(group.test, FALSE) & is.na(group)) {
    warning("group.test will be ignored, since no group is set.")
  }
  if (!identical(group.test, FALSE) & group.long == TRUE) {
    warning("group.test is incompatible with group.long == TRUE and will be ignored.")
  }
  if (!is.logical(factor.numeric) | !is.logical(logical.numeric)) {
    stop("factor.numeric and logical.numeric must each be TRUE or FALSE")
  }
  if (is.logical(group.long)) {
    group.long.align <- "c"
  }
  else if (is.character(group.long)) {
    if (group.long %in% c("l", "r", "c")) {
      group.long.align <- group.long
      group.long <- TRUE
    }
    else {
      stop("group.long must be TRUE, FALSE, or a character l, c, or r.")
    }
  }
  else {
    stop("group.long must be TRUE, FALSE, or a character l, c, or r.")
  }
  if (!is.logical(fixed.digits)) {
    stop("fixed.digits must be TRUE or FALSE.")
  }
  if (fixed.digits) {
    warning("fixed.digits is deprecated and will be removed in a future version in favor of a setting in ")
  }
  if (!is.numeric(col.breaks) & !identical(col.breaks, NA)) {
    stop("col.breaks must be numeric.")
  }
  if (!is.na(group) & !identical(col.breaks, NA) & group.long == 
      FALSE) {
    stop("group cannot be combined with col.breaks unless group.long = TRUE.")
  }
  if (!is.numeric(digits) & !is.list(digits) & !identical(digits, 
                                                          NA)) {
    stop("digits must be numeric.")
  }
  if (!is.list(numformat)) {
    if (length(numformat) > 1) {
      numformat = as.list(numformat)
    }
    else {
      numformat = list(numformat)
    }
  }
  for (fm in 1:length(numformat)) {
    if (is.function(numformat[[fm]])) {
    }
    else if (is.na(numformat[[fm]])) {
      numformat[[fm]] = function(x) x
    }
    else if (is.character(numformat[[fm]])) {
      set_digits <- ifelse(is.na(digits), NULL, digits)
      set_bigmark <- ""
      set_decimalmark <- getOption("OutDec")
      set_percent <- FALSE
      set_prefix = ""
      set_suffix = ""
      if (grepl("comma", numformat[[fm]])) {
        set_bigmark <- ","
        numformat[[fm]] <- gsub("comma", "", numformat[[fm]])
      }
      if (grepl("decimal", numformat[[fm]])) {
        set_bigmark <- "."
        set_decimalmark <- ","
        numformat[[fm]] <- gsub("decimal", "", numformat[[fm]])
      }
      if (grepl("percent", numformat[[fm]])) {
        set_percent <- TRUE
        numformat[[fm]] <- gsub("percent", "", numformat[[fm]])
      }
      if (nchar(numformat[[fm]]) > 0) {
        if (grepl("|", numformat[[fm]])) {
          format_split <- strsplit(numformat[[fm]], 
                                   "|", fixed = TRUE)[[1]]
          set_prefix <- format_split[1]
          set_suffix <- format_split[2]
        }
        else {
          set_prefix <- numformat[[fm]]
        }
      }
      numformat[[fm]] <- formatfunc(percent = set_percent, 
                                    prefix = set_prefix, suffix = set_suffix, digits = digits, 
                                    big.mark = set_bigmark, decimal.mark = set_decimalmark)
    }
    else {
      stop("Each element of numformat must be NA, a string, or a function.")
    }
  }
  if (is.na(obs.function)) {
    obs.function <- "notNA(x)"
    if (!is.null(attr(numformat[[1]], "big.mark"))) {
      obs.function <- paste0("notNA(x, \"", attr(numformat[[1]], 
                                                 "big.mark"), "\")")
    }
    skip.format <- skip.format[!is.na(skip.format)]
    skip.format <- c(skip.format, obs.function)
  }
  if (!is.logical(factor.percent) | !is.logical(factor.counts)) {
    stop("factor.percent and factor.counts must each be TRUE or FALSE.")
  }
  if (!is.character(title)) {
    stop("title must be a character variable.")
  }
  if (!identical(out, NA) & !(out %in% c("viewer", "browser", 
                                         "return", "htmlreturn", "kable", "latex", "latexpage", 
                                         "csv"))) {
    stop("out must be viewer, browser, return, htmlreturn, kable, latex, or latexpage")
  }
  if (identical(out, "csv") & is.na(file)) {
    warning("out = \"csv\" will just return the vtable as a data.frame unless combined with file")
  }
  wts <- NULL
  if (length(group.weights) > 1) {
    wts <- group.weights
  }
  if (length(group.weights) == 1) {
    if (is.character(group.weights)) {
      wts <- data[[group.weights]]
      data[[group.weights]] <- NULL
    }
  }
  if (!identical(group.weights, NA) & is.null(wts)) {
    stop("group.weights must be a vector of length nrow(data), or the name of a column in data")
  }
  if (!is.numeric(wts) & !is.null(wts)) {
    stop("group.weights must be numeric.")
  }
  if ((length(wts) != nrow(data)) & !is.null(wts)) {
    stop("group.weights must be the same length as the number of rows in data")
  }
  if (!is.null(wts)) {
    if (min(wts, na.rm = TRUE) < 0) {
      stop("No negative weights allowed in group.weights")
    }
  }
  if (!is.null(wts)) {
    havewts <- !is.na(wts)
    wts <- wts[havewts]
    data <- subset(data, havewts)
  }
  if (is.matrix(data) & dim(data)[2] == 1) {
    data <- as.data.frame(data)
  }
  var.classes <- sapply(data, function(x) ifelse(is.factor(x), 
                                                 "factor", ifelse(is.logical(x), "logical", ifelse(is.character(x), 
                                                                                                   "character", ifelse(is.numeric(x), "numeric", "other")))))
  labwarning <- FALSE
  for (c in 1:ncol(data)) {
    if (var.classes[c] == "character") {
      if (vtable::nuniq(data[[c]]) <= 6) {
        data[[c]] <- as.factor(data[[c]])
      }
      else {
        if (names(data)[c] %in% vars) {
          warning("You have specified a variable in vars that is a character variable with a large number of different values. It will be excluded. If you are sure you want it in the table, convert it to a factor before calling sumtable.")
        }
        vars <- vars[!(vars == names(data)[c])]
      }
    }
    else if (var.classes[c] == "logical") {
      if (logical.numeric) {
        data[[c]] <- as.numeric(data[[c]])
      }
      else {
        data[[c]] <- factor(data[[c]], levels = c(FALSE, 
                                                  TRUE), labels = logical.labels)
      }
    }
    else if (var.classes[c] == "numeric") {
      if ("labelled" %in% class(data[[c]]) | ("haven_labelled" %in% 
                                              class(data[[c]]) | !is.null(unlist(sjlabelled::get_labels(data[[c]]))))) {
        unlabvals <- length(sjlabelled::get_labels(data[[c]])) == 
          length(sjlabelled::get_labels(data[[c]], non.labelled = TRUE))
        if (!unlabvals) {
          data[[c]] <- as.numeric(data[[c]])
          labwarning <- TRUE
        }
        else {
          suppressWarnings(data[[c]] <- sjlabelled::as_label(data[, 
                                                                  c, drop = FALSE]))
        }
      }
    }
  }
  if (labwarning) {
    warning("Some labelled variables have unlabeled values. Treating these as numeric variables and ignoring labels.")
  }
  var.classes <- sapply(data, function(x) ifelse(is.factor(x), 
                                                 "factor", ifelse(is.logical(x), "logical", ifelse(is.character(x), 
                                                                                                   "character", ifelse(is.numeric(x), "numeric", "other")))))
  factor.warning <- FALSE
  if (any(var.classes == "factor") & !identical(summ, NA) & 
      !factor.numeric) {
    if (is.list(summ) & !identical(col.breaks, NA)) {
      ext.col.breaks <- c(1, col.breaks, ncol(data))
      for (i in 1:length(summ)) {
        if ((!(summ[[i]][1] %in% c("length(x)", obs.function)) | 
             !(summ[[i]][2] %in% "mean(x)")) & any(var.classes[ext.col.breaks[i]:ext.col.breaks[i + 
                                                                                                1]] == "factor")) {
          factor.warning <- TRUE
        }
      }
    }
    else if (!is.list(summ)) {
      if (!(summ[1] %in% c("length(x)", obs.function)) | 
          !(summ[2] %in% "mean(x)")) {
        factor.warning <- TRUE
      }
    }
    else {
      if (!(summ[[1]][1] %in% c("length(x)", obs.function)) | 
          !(summ[[1]][2] %in% "mean(x)")) {
        factor.warning <- TRUE
      }
    }
  }
  if (factor.warning) {
    warning("Factor variables ignore custom summ options. Cols 1 and 2 are count and percentage.\nBeware combining factors with a custom summ unless factor.numeric = TRUE.")
  }
  if (identical(vars, NA)) {
    colkeeps <- sapply(1:ncol(data), function(x) ifelse(is.factor(data[[x]]) | 
                                                          is.numeric(data[[x]]), x, 0))
    if (sum(colkeeps > 0) == 0) {
      stop("It doesn't look like you have any variables that belong in a sumtable. Check your data. Use vars to explicitly choose variables, or convert things to numeric or factor before sending to sumtable.")
    }
    vars <- names(data)[colkeeps[colkeeps > 0]]
    if (!is.na(group)) {
      vars <- vars[vars != group]
    }
    var.classes <- sapply(as.data.frame(data[, vars]), function(x) ifelse(is.factor(x), 
                                                                          "factor", "numeric"))
  }
  else {
    var.classes <- sapply(vars, function(x) ifelse(!(x %in% 
                                                       names(data)), "header", ifelse(is.factor(data[[x]]), 
                                                                                      "factor", ifelse(is.logical(data[[x]]), "logical", 
                                                                                                       ifelse(is.character(data[[x]]), "character", 
                                                                                                              ifelse(is.numeric(data[[x]]), "numeric", "other"))))))
  }
  if (identical(col.breaks, NA)) {
    col.breaks <- length(vars)
  }
  if (utils::tail(col.breaks, 1) < length(vars)) {
    col.breaks[length(col.breaks) + 1] <- length(vars)
  }
  col.windows <- c(0, col.breaks)
  col.vars <- lapply(1:length(col.breaks), function(x) (col.windows[x] + 
                                                          1):col.breaks[x])
  fill.sn <- identical(summ.names, NA)
  if (identical(summ, NA)) {
    summ <- list()
    if (fill.sn) {
      summ.names <- list()
    }
    for (i in 1:length(col.vars)) {
      if (all(var.classes[col.vars[[i]]] == "factor")) {
        summ[[i]] <- c("sum(x)", "mean(x)")
        if (fill.sn & factor.percent) {
          summ.names[[i]] <- c("N", "Percent")
        }
        else {
          summ.names[[i]] <- c("N", "Mean")
        }
        if (!is.null(wts)) {
          summ[[i]] <- c("sum(x)", "stats::weighted.mean(x, w = wts, na.rm = TRUE)")
          if (fill.sn) {
            summ.names[[i]][2] <- paste0(summ.names[[i]][2], 
                                         " (Weighted)")
          }
        }
      }
      else if ((is.na(group) | group.long == TRUE) & length(col.breaks) == 
               1) {
        summ[[i]] <- c(obs.function, "mean(x)", "sd(x)", 
                       "min(x)", "pctile(x)[25]", "pctile(x)[75]", 
                       "max(x)")
        if (fill.sn) {
          summ.names[[i]] <- c("N", "Mean", "Std. Dev.", 
                               "Min", "Pctl. 25", "Pctl. 75", "Max")
        }
        if (add.median) {
          summ[[i]] <- c(obs.function, "mean(x)", "sd(x)", 
                         "min(x)", "pctile(x)[25]", "median(x)", 
                         "pctile(x)[75]", "max(x)")
          if (fill.sn) {
            summ.names[[i]] <- c("N", "Mean", "Std. Dev.", 
                                 "Min", "Pctl. 25", "Pctl. 50", "Pctl. 75", 
                                 "Max")
          }
        }
        if (!is.null(wts)) {
          summ[[i]][summ[[i]] == "mean(x)"] <- "stats::weighted.mean(x, w = wts, na.rm = TRUE)"
          summ[[i]][summ[[i]] == "sd(x)"] <- "weighted.sd(x, w = wts)"
          if (fill.sn) {
            summ.names[[i]][summ.names[[i]] == "Mean"] <- "Wt. Mean"
            summ.names[[i]][summ.names[[i]] == "Std. Dev."] <- "Wt. SD"
          }
        }
      }
      else if ((is.na(group) | group.long == TRUE) & length(col.breaks) > 
               1) {
        summ[[i]] <- c(obs.function, "mean(x)", "sd(x)")
        if (fill.sn) {
          summ.names[[i]] <- c("N", "Mean", "Std. Dev.")
        }
        if (add.median) {
          summ[[i]] <- c(obs.function, "mean(x)", "sd(x)", 
                         "median(x)")
          if (fill.sn) {
            summ.names[[i]] <- c("N", "Mean", "Std. Dev.", 
                                 "Median")
          }
        }
        if (!is.null(wts)) {
          summ[[i]][summ[[i]] == "mean(x)"] <- "stats::weighted.mean(x, w = wts, na.rm = TRUE)"
          summ[[i]][summ[[i]] == "sd(x)"] <- "weighted.sd(x, w = wts)"
          if (fill.sn) {
            summ.names[[i]][summ.names[[i]] == "Mean"] <- "Wt. Mean"
            summ.names[[i]][summ.names[[i]] == "Std. Dev."] <- "Wt. SD"
          }
        }
      }
      else {
        summ[[i]] <- c(obs.function, "mean(x)", "sd(x)")
        if (fill.sn) {
          summ.names[[i]] <- c("N", "Mean", "SD")
        }
        if (add.median) {
          summ[[i]] <- c(obs.function, "mean(x)", "sd(x)", 
                         "median(x)")
          if (fill.sn) {
            summ.names[[i]] <- c("N", "Mean", "SD", 
                                 "Median")
          }
        }
        if (!is.null(wts)) {
          summ[[i]][summ[[i]] == "mean(x)"] <- "stats::weighted.mean(x, w = wts, na.rm = TRUE)"
          summ[[i]][summ[[i]] == "sd(x)"] <- "weighted.sd(x, w = wts)"
          if (fill.sn) {
            summ.names[[i]][summ.names[[i]] == "Mean"] <- "Wt. Mean"
            summ.names[[i]][summ.names[[i]] == "SD"] <- "Wt. SD"
          }
        }
      }
    }
  }
  else if (!is.list(summ)) {
    summ <- lapply(1:length(col.vars), function(x) summ)
  }
  digits.was.list <- is.list(digits)
  if (is.vector(digits)) {
    if (length(digits) > 1) {
      digits.was.list <- TRUE
    }
  }
  if (identical(digits, NA)) {
    digits <- list()
    for (i in 1:length(col.breaks)) {
      digits[[i]] <- rep(3, length(summ[[i]]))
      digits[[i]][1] <- 0
    }
  }
  else if (is.numeric(digits)) {
    if (length(digits) == 1) {
      digopt <- digits
      digits <- list()
      for (i in 1:length(col.breaks)) {
        digits[[i]] <- rep(digopt, length(summ[[i]]))
      }
    }
    else {
      digits <- lapply(1:length(col.breaks), function(x) digits)
    }
  }
  if (fixed.digits & !digits.was.list) {
    for (i in 1:length(summ)) {
      for (j in 1:length(summ[[i]])) {
        calcs <- sapply(vars, function(x) parsefcn_summ(data[[x]], 
                                                        summ[[i]][j]))
        calcs <- calcs[!is.na(calcs)]
        if (is.round(calcs) | summ[[i]][j] == obs.function) {
          digits[[i]][j] <- 0
        }
      }
    }
  }
  if (length(numformat) == 1) {
    single.numformat = numformat[[1]]
    numformat = list()
    for (i in 1:length(vars)) {
      numformat[[i]] = single.numformat
    }
  }
  if (length(numformat) < length(vars)) {
    new.numformat <- list()
    for (v in vars) {
      if (v %in% names(numformat)) {
        new.numformat[[v]] <- numformat[[v]]
      }
      else {
        new.numformat[[v]] <- numformat[[1]]
      }
    }
    numformat <- new.numformat
    rm(new.numformat)
  }
  if (!fill.sn & !is.list(summ.names)) {
    summ.names <- lapply(1:length(col.vars), function(x) summ.names)
  }
  if (identical(summ.names, NA)) {
    summ.names <- list()
    for (i in 1:length(col.vars)) {
      functionsused <- summ[[i]]
      functionsused <- sub("\\(x\\)", "", functionsused)
      firstletters <- toupper(substring(functionsused, 
                                        1, 1))
      summ.names[[i]] <- paste0(firstletters, substring(functionsused, 
                                                        2))
    }
  }
  if (identical(group.test, TRUE)) {
    if (out %in% c("latex", "latexpage") | (isTRUE(getOption("knitr.in.progress")) & 
                                            is.na(out) & isTRUE(knitr::is_latex_output()))) {
      group.test.opts <- list(format = "{name}$={stat}^{{stars}}$")
    }
    else if (out %in% c("return", "kable", "csv") | (isTRUE(getOption("knitr.in.progress")) & 
                                                     is.na(out) & isFALSE(knitr::is_latex_output()) & 
                                                     isFALSE(knitr::is_html_output()))) {
      group.test.opts <- list(format = "{name}={stat}{stars}")
    }
    else {
      group.test.opts <- list(format = "{name}={stat}<sup>{stars}</sup>")
    }
  }
  else if (is.list(group.test)) {
    group.test.opts <- group.test
    group.test <- TRUE
  }
  starnote <- NA_character_
  vartitles <- vars
  grouptitle <- group
  if (identical(labels, TRUE)) {
    labs <- sapply(vars, function(x) attr(data[[x]], "label"))
    has.no.labs <- unlist(sapply(labs, is.null))
    vartitles[!has.no.labs] <- unlist(labs[!has.no.labs])
    if (!is.na(group)) {
      if (!is.null(attr(data[[group]], "label"))) {
        grouptitle <- attr(data[[group]], "label")
      }
    }
  }
  else if (!identical(labels, NA)) {
    if (is.vector(labels)) {
      if (length(labels) == length(vars)) {
        vartitles[!is.na(labels)] <- labels[!is.na(labels)]
      }
      else {
        stop("label vector must have as many elements as there are variables as will be in the sumtable. Use NA elements to fill in, or see help(sumtable) for other label formats that do not require every variable to have a label.")
      }
    }
    else if (dim(labels)[1] > 1 & dim(labels)[2] == 2) {
      temp.df <- data.frame(vars = vars, stringsAsFactors = FALSE)
      labs <- as.data.frame(labels, stringsAsFactors = FALSE)
      names(labs) <- c("vars", "vartitles")
      labs$vars <- as.character(labs$vars)
      labs$vartitles <- as.character(labs$vartitles)
      temp.df$order <- 1:nrow(temp.df)
      temp.df <- merge(temp.df, labs, sort = FALSE, all.x = TRUE)
      temp.df <- temp.df[order(temp.df$order), ]
      temp.df$vartitles[is.na(temp.df$vartitles)] <- temp.df$vars[is.na(temp.df$vartitles)]
      vartitles <- temp.df$vartitles
      if (!is.na(group)) {
        if (sum(labels[[1]] == group) > 0) {
          grouptitle <- labels[labels[[1]] == group, 
                               2]
        }
      }
    }
    else if (dim(labels)[1] == 1 & !is.null(colnames(labels))) {
      labs <- data.frame(vars = colnames(labels), vartitles = as.character(t(labels[1, 
      ])), stringsAsFactors = FALSE)
      temp.df <- data.frame(vars = vars, stringsAsFactors = FALSE)
      temp.df$order <- 1:nrow(temp.df)
      temp.df <- merge(temp.df, labs, sort = FALSE, all.x = TRUE)
      temp.df <- temp.df[order(temp.df$order), ]
      temp.df$vartitles[is.na(temp.df$vartitles)] <- temp.df$vars[is.na(temp.df$vartitles)]
      vartitles <- temp.df$vartitles
      if (!is.na(group)) {
        if (!is.null(labels[[group]])) {
          grouptitle <- labels[[group]][1]
        }
      }
    }
    else {
      stop("Unrecognized label format. See help(vtable).")
    }
  }
  if (is.na(group)) {
    st <- list()
    for (i in 1:length(col.breaks)) {
      st[[i]] <- utils::read.csv(text = paste(c("Variable", 
                                                summ.names[[i]]), collapse = ","), check.names = FALSE)
      contents <- lapply(col.vars[[i]], function(x) {
        summary.row(data, vars[x], st[[i]], vartitles[x], 
                    summ[[i]], var.classes[x], factor.percent, 
                    factor.counts, factor.numeric, digits[[i]], 
                    fixed.digits, wts, numformat[[x]], skip.format, 
                    function(x) eval(parse(text = obs.function)))
      })
      contents <- do.call(rbind, contents)
      st[[i]] <- rbind(st[[i]], contents)
    }
    st <- cbind_unequal(st)
  }
  else if (!group.long) {
    st <- list()
    grouplevels <- sort(unique(data[[group]]))
    for (i in 1:length(grouplevels)) {
      st[[i]] <- utils::read.csv(text = paste(c("Variable", 
                                                summ.names[[1]]), collapse = ","), check.names = FALSE)
      st[[i]][1, ] <- c(paste0("HEADERROW", grouptitle), 
                        paste0(grouplevels[i], "_MULTICOL_c_", length(summ.names[[1]])), 
                        rep("DELETECELL", length(summ.names[[1]]) - 
                              1))
      contents <- lapply(1:length(vars), function(x) summary.row(data[data[[group]] == 
                                                                        grouplevels[i], ], vars[x], st[[i]], vartitles[x], 
                                                                 summ[[1]], var.classes[x], factor.percent, factor.counts, 
                                                                 factor.numeric, digits[[1]], fixed.digits, wts[data[[group]] == 
                                                                                                                  grouplevels[i]], numformat[[x]], skip.format, 
                                                                 function(x) eval(parse(text = obs.function))))
      if (group.test & i == length(grouplevels)) {
        st[[i]] <- utils::read.csv(text = paste(c("Variable", 
                                                  summ.names[[1]], "Test"), collapse = ","), 
                                   check.names = FALSE)
        st[[i]][1, ] <- c(paste0("HEADERROW", grouptitle), 
                          paste0(grouplevels[i], "_MULTICOL_c_", length(summ.names[[1]])), 
                          rep("DELETECELL", length(summ.names[[1]]) - 
                                1), "")
        for (x in 1:length(vars)) {
          test.result <- suppressWarnings(try(independence.test(data[[group]], 
                                                                data[[vars[x]]], w = wts, opts = group.test.opts), 
                                              silent = TRUE))
          if (inherits(test.result, "try-error")) {
            test.result <- ""
          }
          if (!(out %in% c("latex", "latexpage"))) {
            test.result <- gsub("<0", "\\&lt0", test.result)
          }
          contents[[x]]$Test <- c(test.result, rep("", 
                                                   nrow(contents[[x]]) - 1))
        }
      }
      contents <- do.call(rbind, contents)
      st[[i]] <- rbind(st[[i]], contents)
      if (i > 1) {
        st[[i]]$Variable <- NULL
      }
    }
    st <- cbind_unequal(st)
    if (group.test) {
      havenote <- TRUE
      if (!is.null(group.test.opts[["format"]])) {
        havenote <- grepl("\\{stars\\}", group.test.opts[["format"]])
      }
      if (havenote) {
        star.cutoffs <- stat.cutoffs
        star.markers <- c("***", "**", "*")
        if (!is.null(group.test.opts[["star.cutoffs"]])) {
          star.cutoffs <- group.test.opts[["star.cutoffs"]]
        }
        if (!is.null(group.test.opts[["star.markers"]])) {
          star.markers <- group.test.opts[["star.markers"]]
        }
        star.markers <- star.markers[order(-star.cutoffs)]
        star.cutoffs <- star.cutoffs[order(-star.cutoffs)]
        starnote <- paste0(paste0(star.markers, " p<", 
                                  star.cutoffs), collapse = "; ")
        starnote <- paste0("Statistical significance markers: ", 
                           starnote)
      }
    }
  }
  else {
    st <- list()
    grouplevels <- sort(unique(data[[group]]))
    st.all <- list()
    for (j in 1:length(grouplevels)) {
      for (i in 1:length(col.breaks)) {
        st[[i]] <- utils::read.csv(text = paste(c("Variable", 
                                                  summ.names[[i]]), collapse = ","), check.names = FALSE)
        contents <- lapply(col.vars[[i]], function(x) summary.row(data[data[[group]] == 
                                                                         grouplevels[j], ], vars[x], st[[i]], vartitles[x], 
                                                                  summ[[i]], var.classes[x], factor.percent, 
                                                                  factor.counts, factor.numeric, digits[[i]], 
                                                                  fixed.digits, wts[data[[group]] == grouplevels[j]], 
                                                                  numformat[[x]], skip.format, function(x) eval(parse(text = obs.function))))
        summcontents <- do.call(rbind, contents)
        st[[i]] <- rbind(st[[i]], summcontents)
      }
      st.all[[j]] <- cbind_unequal(st)
      header.rows <- st.all[[j]][1, ]
      addrow = 0
      if (j > 1) {
        header.rows[1, ] <- rep("", ncol(header.rows))
        addrow = 1
      }
      header.rows[nrow(header.rows) + addrow, ] <- c(paste0(grouptitle, 
                                                            ": ", grouplevels[j], "_MULTICOL_", group.long.align, 
                                                            "_", ncol(header.rows)), rep("DELETECELL", ncol(header.rows) - 
                                                                                           1))
      st.all[[j]] <- rbind(header.rows, st.all[[j]])
    }
    st <- do.call(rbind, st.all)
  }
  if (identical(col.width, NA) & identical(align, NA)) {
    align <- rep("r", ncol(st))
    align[names(st) == "Variable"] <- "l"
    align[names(st) == "Variable"] <- "@{\\hskip .1in}l"
    if (names(st)[1] == "Variable") {
      align[1] <- "l"
    }
    if (group.test) {
      align[names(st) == "Test"] <- "l"
    }
    align <- paste0(align, collapse = "")
  }
  else {
    align <- paste0("p{", col.width/100, "\\textwidth}")
    if (sum(names(st) == "Variable") > 1) {
      align[names(st) == "Variable"][-1] <- paste0("@{\\hskip .2in}", 
                                                   align[names(st) == "Variable"][-1])
    }
    align <- paste0(align, collapse = "")
  }
  if (identical(col.width, NA)) {
    col.width <- rep(1, ncol(st))
    col.width[names(st) == "Variable"] <- 2
    if (group.test) {
      col.width[names(st) == "Test"] <- 1.5
    }
    totalwidth <- sum(col.width)
    tablescale <- 60 + 20 * (totalwidth >= 2) + 20 * (totalwidth >= 
                                                        3)
    col.width <- (col.width/totalwidth) * tablescale
  }
  if (identical(col.align, NA)) {
    col.align <- rep("right", ncol(st))
    col.align[names(st) == "Variable"] <- "left; padding-left:10px"
    if (names(st)[1] == "Variable") {
      col.align[1] <- "left"
    }
    if (group.test) {
      col.align[names(st) == "Test"] <- "left"
    }
  }
  if (!is.na(group)) {
    names(st)[names(st) != "Variable"] <- paste0(names(st)[names(st) != 
                                                             "Variable"], "_MULTICOL_c_1")
  }
  if (!is.na(note) & !is.na(starnote)) {
    note <- paste0(starnote, ". ", note)
  }
  else if (!is.na(starnote)) {
    note <- starnote
  }
  if (!identical(out, NA) & out %in% c("latex", "latexpage")) {
    if (out == "latex") {
      return(cat(dftoLaTeX(st, file = file, align = align, 
                           anchor = anchor, title = title, note = note, 
                           note.align = note.align, fit.page = fit.page, 
                           no.escape = ifelse(group.test, ncol(st), NA))))
    }
    out.latex <- "\\documentclass{article}\n\\begin{document}\n\n%% sumtable \\{vtable\\}\n\n"
    out.latex <- paste(out.latex, dftoLaTeX(st, align = align, 
                                            anchor = anchor, title = title, note = note, note.align = note.align, 
                                            fit.page = fit.page, no.escape = ifelse(group.test, 
                                                                                    ncol(st), NA)), "\n\n\\end{document}", sep = "")
    if (!is.na(file)) {
      if (!grepl("\\.tex", file)) {
        file <- paste(file, ".tex", sep = "")
      }
      filepath <- file.path(file)
      writeLines(out.latex, filepath)
    }
    return(cat(out.latex))
  }
  out.html <- paste("\n                    <html style=\"font-family:Helvetica,Arial,Sans\">\n                    <head><title>Summary Statistics</title>", 
                    "<style type = \"text/css\">\n                    p {\n                    font-size:smaller;\n                    }\n                    table {\n                    border: 0px;\n                    border-collapse:collapse;\n                    font-size:smaller;\n                    table-layout:fixed;\n                    margin-left:0%;\n                    margin-right:auto;\n                    }\n                    .headtab {\n                    width: 100%;\n                    margin-left:auto;\n                    margin-right:auto;\n                    }\n                    th {\n                    background-color: #FFFFFF;\n                    font-weight:bold;\n                    text-align:left;\n                    }\n                    table tr:nth-child(odd) td {\n                    background-color: #FFFFFF;\n                    padding:4px;\n                    word-wrap: break-word;\n                    word-break:break-all;\n                    }\n                    table tr:nth-child(even) td {\n                    background-color: #D3D3D3;\n                    padding:4px;\n                    word-wrap: break-word;\n                    word-break:break-all;\n                    }</style></head><body>", 
                    sep = "")
  out.html <- paste(out.html, "<table class=\"headtab\">", 
                    "<tr><td style=\"text-align:left\">sumtable {vtable}</td>", 
                    "<td style=\"text-align:right\">Summary Statistics</td></tr></table>", 
                    "<h1>", title, "</h1>")
  out.html <- paste(out.html, dftoHTML(st, out = "htmlreturn", 
                                       col.width = col.width, col.align = col.align, anchor = anchor, 
                                       note = note, note.align = note.align, no.escape = ifelse(group.test, 
                                                                                                ncol(st), NA)), "</body></html>", sep = "")
  if (!is.na(file)) {
    if (identical(out, "csv")) {
      if (!grepl("\\.csv", file)) {
        file <- paste(file, ".csv", sep = "")
      }
      filepath <- file.path(file)
      utils::write.csv(clean_multicol(st), filepath, row.names = FALSE)
    }
    else {
      if (!grepl("\\.htm", file)) {
        file <- paste(file, ".html", sep = "")
      }
      filepath <- file.path(file)
      writeLines(out.html, filepath)
    }
  }
  if (is.na(out)) {
    out = ""
  }
  if (out == "viewer" | out == "browser" | out == "") {
    tempDir <- tempfile()
    dir.create(tempDir)
    htmlpath <- file.path(tempDir, "sumtable.html")
    writeLines(out.html, htmlpath)
  }
  if (out == "kable" | (isTRUE(getOption("knitr.in.progress")) & 
                        out == "")) {
    st <- st[!apply(st, MARGIN = 1, FUN = function(x) !any(!(x == 
                                                               rep("", ncol(st))))), ]
    st <- st[!apply(st, MARGIN = 1, FUN = function(x) propNA(x) == 
                      1), ]
    if (!simple.kable) {
      st <- clean_multicol_kable(st, title, note)
      if (isTRUE(getOption("knitr.in.progress")) & out == 
          "") {
        if (isTRUE(knitr::is_html_output())) {
          st <- kableExtra::kable_styling(st)
        }
      }
      return(st)
    }
    else {
      st <- knitr::kable(clean_multicol(st), caption = title)
      return(st)
    }
  }
  else if (Sys.getenv("RSTUDIO") == "1" & (out == "viewer" | 
                                           out == "")) {
    rstudioapi::viewer(htmlpath)
  }
  else if (Sys.getenv("RSTUDIO") == "" & out == "viewer") {
    stop("out = viewer is not a valid option if RStudio is not running.")
  }
  else if ((Sys.getenv("RSTUDIO") == "" & out == "") | (out == 
                                                        "browser")) {
    utils::browseURL(htmlpath)
  }
  else if (out == "return" | out == "csv") {
    return(clean_multicol(st))
  }
  else if (out == "htmlreturn") {
    return(cat(out.html))
  }
}

summary.row <- function (data, var, st, title, summ, cla, factor.percent, factor.count, 
          factor.numeric, digits, fixed.digits, wts = NULL, fmt = NULL, 
          skip.format = NULL, factor.not.numeric.count = notNA) 
{
  numcols <- length(summ)
  if (cla == "header") {
    st[1, ] <- c(paste0(title, "_MULTICOL_l_all"), rep("DELETECELL", 
                                                       numcols))
  }
  else if (cla == "factor" & !factor.numeric) {
    va <- data[[var]]
    nonmiss <- factor.not.numeric.count(va)
    nonmissnum <- notNA(va)
    st[1, ] <- c(title, nonmiss, "", rep("", numcols - 2))
    mat <- as.data.frame(table(va))
    if (nonmissnum > 0) {
      matlabel <- stats::aggregate(y ~ x, data.frame(y = 1, 
                                                     x = va), FUN = factor.not.numeric.count, drop = FALSE)
      matlabel$y[is.na(matlabel$y)] <- 0
      propcalc <- mat$Freq/nonmissnum
      if (!is.null(wts) & grepl("wts", summ[2])) {
        propcalc <- sapply(mat$va, function(x) stats::weighted.mean(va == 
                                                                      x, w = wts, na.rm = TRUE))
      }
    }
    else {
      matlabel <- stats::aggregate(y ~ x, data.frame(y = 1, 
                                                     x = factor(levels(va), levels = levels(va))), 
                                   FUN = factor.not.numeric.count, drop = FALSE)
      matlabel$y <- 0
      propcalc <- rep(NA, length(mat$Freq))
    }
    propcalc <- propcalc * (100^factor.percent)
    mat$va <- paste("...", mat$va)
    mat$Freq <- matlabel$y
    if (fixed.digits) {
      mat$Prop <- sapply(1:length(propcalc), function(x) format(propcalc[x], 
                                                                digits = max(digits[2] - 2 * factor.percent, 
                                                                             1), nsmall = max(digits[2] - 2 * factor.percent, 
                                                                                              0), scientific = FALSE))
      st[1, 2] <- format(as.numeric(st[1, 2]), max(digits = digits[1], 
                                                   1), nsmall = digits[1], scientific = FALSE)
    }
    else {
      mat$Prop <- round(propcalc, digits = max(digits[2] - 
                                                 2 * factor.percent, 0))
    }
    if (!factor.count) {
      mat$Freq <- ""
    }
    if (factor.percent) {
      mat$Prop <- paste0(mat$Prop, "%")
    }
    if (ncol(mat) < ncol(st)) {
      mat[, (ncol(mat) + 1):(ncol(st))] <- ""
    }
    if (nonmissnum == 0) {
      mat$Prop <- ""
    }
    names(mat) <- names(st)
    st <- rbind(st, mat)
  }
  else if (cla == "factor" & factor.numeric) {
    st[1, ] <- c(title, rep("", numcols))
    va <- data[, var]
    mat <- stats::model.matrix(~-1 + va)
    facnames <- paste("...", levels(va))
    results <- lapply(1:ncol(mat), function(x) lapply(summ, 
                                                      function(y) parsefcn_summ(mat[, x], y, wts = wts[!is.na(va)])))
    results <- lapply(1:length(results), function(x) sapply(1:length(results[[x]]), 
                                                            function(y) ifelse(summ[y] %in% skip.format, results[[x]][[y]], 
                                                                               fmt(results[[x]][[y]]))))
    if (paste0(deparse(fmt), collapse = "") == "function (x) x") {
      if (fixed.digits) {
        results <- lapply(results, function(x) sapply(1:length(x), 
                                                      function(y) ifelse(is.character(x[y]), x[y], 
                                                                         format(x[y], digits = digits[y], nsmall = max(digits[y], 
                                                                                                                       1), scientific = FALSE))))
      }
      else {
        results <- lapply(results, function(x) sapply(1:length(x), 
                                                      function(y) ifelse(is.character(x[y]), x[y], 
                                                                         as.character(round(x[y], digits = digits[y])))))
      }
    }
    results <- lapply(1:length(results), function(x) c(facnames[x], 
                                                       results[[x]]))
    results <- as.data.frame(do.call(rbind, results))
    names(results) <- names(st)
    st <- rbind(st, results)
  }
  else {
    va <- data[[var]]
    results <- lapply(summ, function(y) parsefcn_summ(va, 
                                                      y, wts = wts[!is.na(va)]))
    results <- lapply(1:length(results), function(y) ifelse(summ[y] %in% 
                                                              skip.format, results[[y]], fmt(results[[y]])))
    if (paste0(deparse(fmt), collapse = "") == "function (x) x") {
      if (fixed.digits) {
        results <- sapply(1:length(results), function(y) ifelse(is.character(results[[y]]), 
                                                                results[[y]], format(results[[y]], digits = max(digits[y], 
                                                                                                                1), nsmall = digits[y], scientific = FALSE)))
      }
      else {
        results <- sapply(1:length(results), function(y) ifelse(is.character(results[[y]]), 
                                                                results[[y]], round(results[[y]], digits = max(digits[y], 
                                                                                                               1))))
      }
    }
    st[1, ] <- c(title, results)
  }
  return(st)
}

parsefcn_summ <- function (x, y, ...) 
{
  list2env(list(...), envir = environment())
  if (!any(sapply(c("anyNA", "propNA", "countNA"), function(z) grepl(z, 
                                                                     y)))) {
    x <- x[!is.na(x)]
  }
  result <- suppressWarnings(try(eval(parse(text = y)), silent = TRUE))
  if (inherits(result, "try-error")) {
    result <- NA
  }
  return(result)
}

cbind_unequal <- function (x) 
{
  rowseach <- sapply(x, nrow)
  mostrows <- max(rowseach)
  for (i in 1:length(x)) {
    if (rowseach[i] < mostrows) {
      x[[i]][(rowseach[i] + 1):mostrows, ] <- ""
    }
  }
  return(do.call(cbind, x))
}
