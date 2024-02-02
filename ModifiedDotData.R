library(tidyr)

PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
}


# ModifiedDotData <- function(
#     object,
#     group.by,
#     markers,
#     cols.use = c("lightgrey", "blue"),
#     col.min = -2.5,
#     col.max = 2.5,
#     dot.min = 0,
#     dot.scale = 6,
#     scale.by = 'radius'
# ) {
#     assay <- assay %||% DefaultAssay(object = object)
#     DefaultAssay(object = object) <- assay
#     scale.func <- switch(
#         EXPR = scale.by,
#         'size' = scale_size,
#         'radius' = scale_radius,
#         stop("'scale.by' must be either 'size' or 'radius'")
#     )
#     if (!missing(x = group.by)) {
#         object <- SetAllIdent(object = object, id = group.by)
#     }
#     
#     genes.plot <- intersect(row.names(object@raw.data), markers$Markers)
#     
#     data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
#     colnames(x = data.to.plot) <- genes.plot
#     data.to.plot$cell <- rownames(x = data.to.plot)
#     data.to.plot$id <- object@ident
#     
#     data.to.plot %>% gather(
#         key = genes.plot,
#         value = expression,
#         -c(cell, id)
#     ) -> data.to.plot
#     data.to.plot %>%
#         group_by(id, genes.plot) %>%
#         summarize(
#             avg.exp = mean(expm1(x = expression)),
#             pct.exp = PercentAbove(x = expression, threshold = 0)
#         ) -> data.to.plot
#     data.to.plot %>%
#         ungroup() %>%
#         group_by(genes.plot) %>%
#         mutate(avg.exp.scale = scale(x = avg.exp)) %>%
#         mutate(avg.exp.scale = MinMax(
#             data = avg.exp.scale,
#             max = col.max,
#             min = col.min
#         )) ->  data.to.plot
#     data.to.plot$genes.plot <- factor(
#         x = data.to.plot$genes.plot,
#         levels = rev(x = genes.plot)
#     )
#     # data.to.plot$genes.plot <- factor(
#     #   x = data.to.plot$genes.plot,
#     #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
#     # )
#     data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
#     data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
#     data.to.plot <- merge(data.to.plot, markers, by.x = "genes.plot", by.y = "Markers")
#     
#     return(data.to.plot)
# }


ModifyDotData <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    group.by = NULL,
    split.by = NULL,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
) {
    if (is.null(assay)) {
        assay = DefaultAssay(object = object)
    }
    
    features = intersect(features, rownames(object))
    
    DefaultAssay(object = object) <- assay
    scale.func <- switch(
        EXPR = scale.by,
        'size' = scale_size,
        'radius' = scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )
    data.features <- FetchData(object = object, vars = features)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    } else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = '_')
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(
        X = unique(x = data.features$id),
        FUN = function(ident) {
            data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
            avg.exp <- apply(
                X = data.use,
                MARGIN = 2,
                FUN = function(x) {
                    return(mean(x = expm1(x = x)))
                }
            )
            pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
            return(list(avg.exp = avg.exp, pct.exp = pct.exp))
        }
    )
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(
        X = names(x = data.plot),
        FUN = function(x) {
            data.use <- as.data.frame(x = data.plot[[x]])
            data.use$features.plot <- rownames(x = data.use)
            data.use$id <- x
            return(data.use)
        }
    )
    data.plot <- do.call(what = 'rbind', args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(
        X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min, max = col.max)
            } else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        }
    )
    
    
    
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(
        x = data.plot$features.plot,
        levels = rev(x = features)
    )
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(
            X = as.character(x = data.plot$id),
            FUN = gsub,
            FUN.VALUE = character(length = 1L),
            pattern =  paste0(
                '^((',
                paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
                ')_)'
            ),
            replacement = '',
            USE.NAMES = FALSE
        )
        data.plot$colors <- mapply(
            FUN = function(color, value) {
                return(colorRampPalette(colors = c('grey', color))(20)[value])
            },
            color = cols[splits.use],
            value = avg.exp.scaled
        )
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = 'avg.exp.scaled', no = 'colors')
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
    }
    
    return(data.plot)
}

