trellis.par.set(canonical.theme("pdf"))

font.settings <- list(font = 1)

ps.theme.shared <- list(strip.background = list(alpha = 1, col = "white"),
                        layout.heights = list(strip = 1.5),
                        par.xlab.text = font.settings
                        )

ps.cex = .8

ps.theme.1 <- c(ps.theme.shared,
                list(superpose.symbol = list(
                       cex = ps.cex,
                       pch = c(21, 21),
                       col = c( "black", "black"),
                       fill = c( "white", "black"))
                     )
                )


ps.theme.1 <- c(ps.theme.shared,
                list(superpose.symbol = list(
                       cex = ps.cex,
                       pch = c(19, 22),
                       col = c( "#0080ff", "black"),
                       fill = c( "#0080ff", "black")),
                     superpose.line = list(
                       col = c( "#0080ff", "black")
                       )
                     )
                )

model.cat.col <- c("exp" = "#0080ff", "naive" = "#8a8a8a",
                   "basic" = '#E5AA00', "M" = '#E5AA00',
                   "s" = '#E7298A', "Ms" = '#E7298A')
model.cat.symbol <- c("exp" = 19, "naive" = 19,
                      "basic" = 19, "M" = 15,
                      "s" = 19, "Ms" = 15)
model.cat.lty <- c("exp" = 1, "naive" = 1,
                      "basic" = 1, "M" = 1,
                      "s" = 2, "Ms" = 2)
ps.theme.2 <- c(ps.theme.shared,
                list(superpose.symbol = list(
                       cex = ps.cex,
                       pch = model.cat.symbol[levels(pssbo.ci$model.cat)],
                       col = model.cat.col[levels(pssbo.ci$model.cat)]#,
                       #fill = model.cat.col[levels(pssbo.ci$model.cat)]
                       ),
                     superpose.line = list(
                       lty = model.cat.lty,
                       col = model.cat.col[levels(pssbo.ci$model.cat)]
                       )
                     )
                )

ps.theme.2.m <- c(ps.theme.shared,
                list(superpose.symbol = list(
                       cex = ps.cex,
                       col = model.cat.col[c(5, 5, 5, 5, 3, 3, 3, 3, 2, 1)]
                       ),
                     superpose.line = list(
                       cex = ps.cex,
                       col = model.cat.col[c(5, 5, 5, 5, 3, 3, 3, 3, 2, 1)]
                       )
                     )
                )

ps.theme.2.m.noexp <- c(ps.theme.shared,
                list(superpose.symbol = list(
                       cex = ps.cex,
                       col = model.cat.col[c(5, 5, 5, 5, 3, 3, 3, 3)]
                       ),
                     superpose.line = list(
                       cex = ps.cex,
                       col = model.cat.col[c(5, 5, 5, 5, 3, 3, 3, 3)]
                       )
                     )
                )

ps.theme.3 <- c(ps.theme.shared,
                list(superpose.symbol = list(pch = c(19),
                       col = rev(c("black", "red")),
                       fill = rev(c("black", "red")))
                     )
                )


exposed.col <- "#067621"
ps.theme.ps <- c(
  ps.theme.shared,
  list(superpose.line = list(
         col = c("black", exposed.col)
         ),
       superpose.symbol = list(
         col = c("black", exposed.col),
         pch = c(1, 4),
         cex = .6
         )
       )
  )


#trellis.par.set(theme = ps.theme.shared)

#col = c("#0080ff", '#E7298A', '#E6AB02')
#fill = c("#0080ff", '#E7298A', '#E6AB02')


# labels

p.0.expr <- expression(hat(italic(p)) ^(0))
atet.expr <- expression(hat(delta))
atet.error.expr <- expression(paste("Error in ", hat(delta)))
atet.rel.error.expr <- expression(paste("Relative error in ", hat(delta)))
rel.risk.expr <- expression("Estimated relative risk")
rel.risk.rel.error.expr <- expression("Percent error in relative risk")
or.expr <- expression("Estimated odds ratio")
percent.poss.expr <- expression(paste("Percent of possible error in ", hat(delta)))
percent.red.expr <- expression(paste("Percent error reduction for ", delta))
percent.red.sds.expr <- expression(paste("Percent remaining error reduction for ", delta))
