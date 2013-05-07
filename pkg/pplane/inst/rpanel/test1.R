#install.packages("rpanel")
library(rpanel)

# textbox
plotf <- function(panel) {
    with(panel, {
                pars   <- as.numeric(pars)
                xgrid <- seq(0.1, max(c(pars[3], 5), na.rm = TRUE), length = 50)
                dgrid <- df(xgrid, pars[1], pars[2])
                plot(xgrid, dgrid, type = "l", col = "blue", lwd = 3)
                if (!is.na(pars[3])) {
                    lines(rep(pars[3], 2), c(0, 0.95 * max(dgrid)), lty = 2, col = "red")
                    text(pars[3], max(dgrid), as.character(pars[3]), col = "red")
                }
            })
    panel
}

panel <- rp.control(pars = c(5, 10, NA))
rp.textentry(panel, pars, plotf, labels = c("df1", "df2", "observed"),
        initval = c(10, 5, 3))
rp.do(panel, plotf)


# listbox
data.plotfn <- function(panel) {
    if (panel$plot.type == "histogram"){
        hist(panel$x)
    }else if (panel$plot.type == "boxplot"){
        boxplot(panel$x)
    }else
        plot(density(panel$x))
    panel
}
panel <- rp.control(x = rnorm(50))
rp.listbox(panel, plot.type,
        c("histogram", "boxplot", "density estimate"),
        action = data.plotfn, title = "Plot type", name="LBPlotType") 

rp.widget.dispose(panel,LBPlotType)

# can now add with different entries
rp.listbox(panel, plot.type,
        c("histogram", "density estimate", "boxplot"),
        action = data.plotfn, title = "Plot type", name="LBPlotType") 
