# R-Code-For-Milk-Lab
R version 4.6.1 (2026-06-24) -- "Happy Hop"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin23

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/Lab 2 Calibration Curve/.RData]

> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # Calibration curves for all 6 color channels:
> #   Red, Green, Blue, Hue, Saturation, Lightness
> #
> # Data: ONE Color Grab reading per protein standard (original readings).
> # Standards: 1.0, 2.0, 3.0, 4.0 mg/mL dilutions + 5.0 mg/mL stock
> #            (80:20 casein:BSA, from the lab procedure).
> #
> # Each graph has: axis labels, title, linear regression line,
> # equation, R^2, and error bars.
> #
> # ERROR BARS: each standard was read only once, so there is no
> # replicate standard deviation. The error bars show +/- the standard
> # error of the regression (s_y/x) for that channel, which is the typical
> # scatter of a single reading around the calibration line. Say this
> # in your figure caption.
> #
> # HOW TO USE: open in R, click "Source". Files are saved to the folder
> # printed at the end.
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # ---- 2. HEX CODE -> R, G, B, HUE, SATURATION, LIGHTNESS ---------------
> # R, G, B: 0-255.  Hue: degrees (0-360).  Saturation, Lightness: percent.
> rgb_vals <- t(col2rgb(standards$hex))
> r <- rgb_vals[, 1] / 255
> g <- rgb_vals[, 2] / 255
> b <- rgb_vals[, 3] / 255
> cmax  <- pmax(r, g, b)
> cmin  <- pmin(r, g, b)
> delta <- cmax - cmin
> 
> L <- (cmax + cmin) / 2
> S <- ifelse(delta == 0, 0, delta / (1 - abs(2 * L - 1)))
> H <- ifelse(delta == 0, 0,
+             ifelse(cmax == r, 60 * (((g - b) / delta) %% 6),
+                    ifelse(cmax == g, 60 * ((b - r) / delta + 2),
+                           60 * ((r - g) / delta + 4))))
> 
> standards$Red        <- rgb_vals[, 1]
> standards$Green      <- rgb_vals[, 2]
> standards$Blue       <- rgb_vals[, 3]
> standards$Hue        <- H
> standards$Saturation <- S * 100
> standards$Lightness  <- L * 100
> 
> cat("\n===== Color values for each standard =====\n")

===== Color values for each standard =====
> print(standards, row.names = FALSE, digits = 4)
 conc_mg_mL     hex Red Green Blue   Hue Saturation Lightness
          1 #EABB44 234   187   68 43.01      79.81     59.22
          2 #E2AD3D 226   173   61 40.73      73.99     56.27
          3 #DF9E37 223   158   55 36.79      72.41     54.51
          4 #B28E38 178   142   56 42.30      52.14     45.88
          5 #AB9321 171   147   33 49.57      67.65     40.00
> write.csv(standards, "standards_color_values.csv", row.names = FALSE)
> 
> # ---- 3. AXIS LABELS FOR EACH CHANNEL ----------------------------------
> channels <- c("Red", "Green", "Blue", "Hue", "Saturation", "Lightness")
> y_labels <- c(Red        = "Red channel intensity (0-255)",
+               Green      = "Green channel intensity (0-255)",
+               Blue       = "Blue channel intensity (0-255)",
+               Hue        = "Hue (degrees)",
+               Saturation = "Saturation (%)",
+               Lightness  = "Lightness (%)")
> 
> # ---- 4. ONE GRAPH PER CHANNEL -----------------------------------------
> results <- data.frame()
> 
> for (i in seq_along(channels)) {
+     ch  <- channels[i]
+     dat <- data.frame(conc = standards$conc_mg_mL,
+                       y    = standards[[ch]],
+                       hex  = standards$hex)
+     
+     # Linear regression: channel value vs. protein concentration
+     fit   <- lm(y ~ conc, data = dat)
+     s     <- summary(fit)
+     slope <- coef(fit)[2]
+     int   <- coef(fit)[1]
+     r2    <- s$r.squared
+     s_yx  <- s$sigma                      # standard error of the regression
+     
+     results <- rbind(results, data.frame(
+         channel       = ch,
+         slope         = slope,
+         slope_SE      = s$coefficients[2, 2],
+         intercept     = int,
+         intercept_SE  = s$coefficients[1, 2],
+         R_squared     = r2,
+         s_yx          = s_yx,
+         row.names     = NULL))
+     
+     # Equation label, e.g.  y = -4.882x + 65.82,  R^2 = 0.936
+     eq <- sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+                   slope, ifelse(int < 0, "-", "+"), abs(int), r2)
+     
+     # Put the equation in an empty corner (bottom-left if the line goes
+     # down, top-left if it goes up)
+     if (slope < 0) {
+         eq_y <- min(dat$y - s_yx); eq_vjust <- 0
+     } else {
+         eq_y <- max(dat$y + s_yx); eq_vjust <- 1
+     }
+     
+     p <- ggplot(dat, aes(x = conc, y = y)) +
+         geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
+                     colour = "black", linewidth = 0.5) +
+         geom_errorbar(aes(ymin = y - s_yx, ymax = y + s_yx),
+                       width = 0.15, linewidth = 0.4) +
+         # each point is filled with the actual color of that disk
+         geom_point(aes(fill = hex), shape = 21, size = 2.6,
+                    colour = "black", stroke = 0.6) +
+         scale_fill_identity() +
+         annotate("text", x = min(dat$conc), y = eq_y, hjust = 0,
+                  vjust = eq_vjust, size = 2.6, parse = TRUE, label = eq) +
+         scale_x_continuous(breaks = standards$conc_mg_mL) +
+         labs(title = paste0("Calibration curve: ", tolower(ch),
+                             " vs. protein concentration"),
+              x = "Protein standard concentration (mg/mL)",
+              y = y_labels[ch]) +
+         theme_classic(base_size = 8) +
+         theme(plot.title = element_text(size = 8, face = "bold"),
+               axis.text  = element_text(colour = "black"))
+     
+     print(p)
+     file_name <- sprintf("Figure_%d_%s_calibration.png", i, ch)
+     ggsave(file_name, p, width = 3.3, height = 2.8, units = "in", dpi = 600)
+     cat("Saved", file_name, "\n")
+ }
Saved Figure_1_Red_calibration.png 
Saved Figure_2_Green_calibration.png 
Saved Figure_3_Blue_calibration.png 
Saved Figure_4_Hue_calibration.png 
Saved Figure_5_Saturation_calibration.png 
Saved Figure_6_Lightness_calibration.png 
> 
> # ---- 5. SUMMARY TABLE: WHICH CHANNEL IS MOST LINEAR? -------------------
> results <- results[order(-results$R_squared), ]
> cat("\n===== Linear fits for all 6 channels (best R^2 first) =====\n")

===== Linear fits for all 6 channels (best R^2 first) =====
> print(results, row.names = FALSE, digits = 4)
    channel   slope slope_SE intercept intercept_SE R_squared   s_yx
  Lightness  -4.882   0.7371     65.82        2.445    0.9360  2.331
      Green -11.100   2.2591    194.70        7.492    0.8895  7.144
        Red -17.400   3.8884    258.60       12.897    0.8697 12.296
       Blue  -7.500   2.0551     77.10        6.816    0.8162  6.499
 Saturation  -4.618   2.7456     83.05        9.106    0.4853  8.682
        Hue   1.467   1.4661     38.07        4.862    0.2503  4.636
> write.csv(results, "calibration_fits_all_channels.csv", row.names = FALSE)
> 
> cat("\nBest channel (highest R^2):", results$channel[1],
+     sprintf("(R^2 = %.3f)", results$R_squared[1]), "\n")

Best channel (highest R^2): Lightness (R^2 = 0.936) 
> cat("\nAll files saved to:", getwd(), "\n")

All files saved to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # Calibration curves for all 6 color channels:
> #   Red, Green, Blue, Hue, Saturation, Lightness
> #
> # Data: ONE Color Grab reading per protein standard (original readings).
> # Standards: 1.0, 2.0, 3.0, 4.0 mg/mL dilutions + 5.0 mg/mL stock
> #            (80:20 casein:BSA, from the lab procedure).
> #
> # Each graph has: axis labels, title, linear regression line,
> # equation, and R^2. No error bars. Each graph is drawn in its own color.
> #
> # HOW TO USE: open in R, click "Source". Files are saved to the folder
> # printed at the end.
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # ---- 2. GRAPH COLORS (change any of these if you want) ------------------
> graph_colors <- c(Red        = "#E34948",   # red
+                   Green      = "#008300",   # green
+                   Blue       = "#2A78D6",   # blue
+                   Hue        = "#6250D6",   # violet
+                   Saturation = "#EB6834",   # orange
+                   Lightness  = "#1BAF7A")   # teal
> 
> # ---- 3. HEX CODE -> R, G, B, HUE, SATURATION, LIGHTNESS ---------------
> # R, G, B: 0-255.  Hue: degrees (0-360).  Saturation, Lightness: percent.
> rgb_vals <- t(col2rgb(standards$hex))
> r <- rgb_vals[, 1] / 255
> g <- rgb_vals[, 2] / 255
> b <- rgb_vals[, 3] / 255
> cmax  <- pmax(r, g, b)
> cmin  <- pmin(r, g, b)
> delta <- cmax - cmin
> 
> L <- (cmax + cmin) / 2
> S <- ifelse(delta == 0, 0, delta / (1 - abs(2 * L - 1)))
> H <- ifelse(delta == 0, 0,
+             ifelse(cmax == r, 60 * (((g - b) / delta) %% 6),
+                    ifelse(cmax == g, 60 * ((b - r) / delta + 2),
+                           60 * ((r - g) / delta + 4))))
> 
> standards$Red        <- rgb_vals[, 1]
> standards$Green      <- rgb_vals[, 2]
> standards$Blue       <- rgb_vals[, 3]
> standards$Hue        <- H
> standards$Saturation <- S * 100
> standards$Lightness  <- L * 100
> 
> cat("\n===== Color values for each standard =====\n")

===== Color values for each standard =====
> print(standards, row.names = FALSE, digits = 4)
 conc_mg_mL     hex Red Green Blue   Hue Saturation Lightness
          1 #EABB44 234   187   68 43.01      79.81     59.22
          2 #E2AD3D 226   173   61 40.73      73.99     56.27
          3 #DF9E37 223   158   55 36.79      72.41     54.51
          4 #B28E38 178   142   56 42.30      52.14     45.88
          5 #AB9321 171   147   33 49.57      67.65     40.00
> write.csv(standards, "standards_color_values.csv", row.names = FALSE)
> 
> # ---- 4. AXIS LABELS FOR EACH CHANNEL ----------------------------------
> channels <- c("Red", "Green", "Blue", "Hue", "Saturation", "Lightness")
> y_labels <- c(Red        = "Red channel intensity (0-255)",
+               Green      = "Green channel intensity (0-255)",
+               Blue       = "Blue channel intensity (0-255)",
+               Hue        = "Hue (degrees)",
+               Saturation = "Saturation (%)",
+               Lightness  = "Lightness (%)")
> 
> # ---- 5. ONE GRAPH PER CHANNEL -----------------------------------------
> results <- data.frame()
> 
> for (i in seq_along(channels)) {
+     ch  <- channels[i]
+     col <- graph_colors[ch]
+     dat <- data.frame(conc = standards$conc_mg_mL, y = standards[[ch]])
+     
+     # Linear regression: channel value vs. protein concentration
+     fit   <- lm(y ~ conc, data = dat)
+     s     <- summary(fit)
+     slope <- coef(fit)[2]
+     int   <- coef(fit)[1]
+     r2    <- s$r.squared
+     
+     results <- rbind(results, data.frame(
+         channel       = ch,
+         slope         = slope,
+         slope_SE      = s$coefficients[2, 2],
+         intercept     = int,
+         intercept_SE  = s$coefficients[1, 2],
+         R_squared     = r2,
+         regression_SE = s$sigma,
+         row.names     = NULL))
+     
+     # Equation label, e.g.  y = -4.882x + 65.82,  R^2 = 0.936
+     eq <- sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+                   slope, ifelse(int < 0, "-", "+"), abs(int), r2)
+     
+     # Put the equation in an empty corner (bottom-left if the line goes
+     # down, top-left if it goes up)
+     y_pad <- 0.08 * diff(range(dat$y))
+     if (slope < 0) {
+         eq_y <- min(dat$y) - y_pad; eq_vjust <- 0
+     } else {
+         eq_y <- max(dat$y) + y_pad; eq_vjust <- 1
+     }
+     
+     p <- ggplot(dat, aes(x = conc, y = y)) +
+         geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
+                     colour = col, linewidth = 0.7) +
+         geom_point(shape = 21, size = 2.8, fill = col,
+                    colour = "white", stroke = 0.8) +
+         annotate("text", x = min(dat$conc), y = eq_y, hjust = 0,
+                  vjust = eq_vjust, size = 2.6, parse = TRUE, label = eq) +
+         scale_x_continuous(breaks = standards$conc_mg_mL) +
+         labs(title = paste0("Calibration curve: ", tolower(ch),
+                             " vs. protein concentration"),
+              x = "Protein standard concentration (mg/mL)",
+              y = y_labels[ch]) +
+         theme_classic(base_size = 8) +
+         theme(plot.title = element_text(size = 8, face = "bold"),
+               axis.text  = element_text(colour = "black"))
+     
+     print(p)
+     file_name <- sprintf("Figure_%d_%s_calibration.png", i, ch)
+     ggsave(file_name, p, width = 3.3, height = 2.8, units = "in", dpi = 600)
+     cat("Saved", file_name, "\n")
+ }
Saved Figure_1_Red_calibration.png 
Saved Figure_2_Green_calibration.png 
Saved Figure_3_Blue_calibration.png 
Saved Figure_4_Hue_calibration.png 
Saved Figure_5_Saturation_calibration.png 
Saved Figure_6_Lightness_calibration.png 
> 
> # ---- 6. SUMMARY TABLE: WHICH CHANNEL IS MOST LINEAR? -------------------
> results <- results[order(-results$R_squared), ]
> cat("\n===== Linear fits for all 6 channels (best R^2 first) =====\n")

===== Linear fits for all 6 channels (best R^2 first) =====
> print(results, row.names = FALSE, digits = 4)
    channel   slope slope_SE intercept intercept_SE R_squared regression_SE
  Lightness  -4.882   0.7371     65.82        2.445    0.9360         2.331
      Green -11.100   2.2591    194.70        7.492    0.8895         7.144
        Red -17.400   3.8884    258.60       12.897    0.8697        12.296
       Blue  -7.500   2.0551     77.10        6.816    0.8162         6.499
 Saturation  -4.618   2.7456     83.05        9.106    0.4853         8.682
        Hue   1.467   1.4661     38.07        4.862    0.2503         4.636
> write.csv(results, "calibration_fits_all_channels.csv", row.names = FALSE)
> 
> cat("\nBest channel (highest R^2):", results$channel[1],
+     sprintf("(R^2 = %.3f)", results$R_squared[1]), "\n")

Best channel (highest R^2): Lightness (R^2 = 0.936) 
> cat("\nAll files saved to:", getwd(), "\n")

All files saved to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # Calibration curves for all 6 color channels:
> #   Red, Green, Blue, Hue, Saturation, Lightness
> #
> # Data: ONE Color Grab reading per protein standard (original readings).
> # Standards: 1.0, 2.0, 3.0, 4.0 mg/mL dilutions + 5.0 mg/mL stock
> #            (80:20 casein:BSA, from the lab procedure).
> #
> # Each graph has: axis labels, title, linear regression line,
> # equation, and R^2. No error bars.
> # Colors: every line of best fit is black. Dots are red, green, and blue
> # on the Red, Green, and Blue graphs; dots are black on the Hue,
> # Saturation, and Lightness graphs.
> #
> # HOW TO USE: open in R, click "Source". Files are saved to the folder
> # printed at the end.
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # ---- 2. GRAPH COLORS (change any of these if you want) ------------------
> # Dot color for each graph:
> dot_colors <- c(Red        = "#E34948",   # red
+                 Green      = "#008300",   # green
+                 Blue       = "#2A78D6",   # blue
+                 Hue        = "black",
+                 Saturation = "black",
+                 Lightness  = "black")
> 
> # Line of best fit color (same on every graph):
> line_color <- "black"
> 
> # ---- 3. HEX CODE -> R, G, B, HUE, SATURATION, LIGHTNESS ---------------
> # R, G, B: 0-255.  Hue: degrees (0-360).  Saturation, Lightness: percent.
> rgb_vals <- t(col2rgb(standards$hex))
> r <- rgb_vals[, 1] / 255
> g <- rgb_vals[, 2] / 255
> b <- rgb_vals[, 3] / 255
> cmax  <- pmax(r, g, b)
> cmin  <- pmin(r, g, b)
> delta <- cmax - cmin
> 
> L <- (cmax + cmin) / 2
> S <- ifelse(delta == 0, 0, delta / (1 - abs(2 * L - 1)))
> H <- ifelse(delta == 0, 0,
+             ifelse(cmax == r, 60 * (((g - b) / delta) %% 6),
+                    ifelse(cmax == g, 60 * ((b - r) / delta + 2),
+                           60 * ((r - g) / delta + 4))))
> 
> standards$Red        <- rgb_vals[, 1]
> standards$Green      <- rgb_vals[, 2]
> standards$Blue       <- rgb_vals[, 3]
> standards$Hue        <- H
> standards$Saturation <- S * 100
> standards$Lightness  <- L * 100
> 
> cat("\n===== Color values for each standard =====\n")

===== Color values for each standard =====
> print(standards, row.names = FALSE, digits = 4)
 conc_mg_mL     hex Red Green Blue   Hue Saturation Lightness
          1 #EABB44 234   187   68 43.01      79.81     59.22
          2 #E2AD3D 226   173   61 40.73      73.99     56.27
          3 #DF9E37 223   158   55 36.79      72.41     54.51
          4 #B28E38 178   142   56 42.30      52.14     45.88
          5 #AB9321 171   147   33 49.57      67.65     40.00
> write.csv(standards, "standards_color_values.csv", row.names = FALSE)
> 
> # ---- 4. AXIS LABELS FOR EACH CHANNEL ----------------------------------
> channels <- c("Red", "Green", "Blue", "Hue", "Saturation", "Lightness")
> y_labels <- c(Red        = "Red channel intensity (0-255)",
+               Green      = "Green channel intensity (0-255)",
+               Blue       = "Blue channel intensity (0-255)",
+               Hue        = "Hue (degrees)",
+               Saturation = "Saturation (%)",
+               Lightness  = "Lightness (%)")
> 
> # ---- 5. ONE GRAPH PER CHANNEL -----------------------------------------
> results <- data.frame()
> 
> for (i in seq_along(channels)) {
+     ch  <- channels[i]
+     col <- dot_colors[ch]
+     dat <- data.frame(conc = standards$conc_mg_mL, y = standards[[ch]])
+     
+     # Linear regression: channel value vs. protein concentration
+     fit   <- lm(y ~ conc, data = dat)
+     s     <- summary(fit)
+     slope <- coef(fit)[2]
+     int   <- coef(fit)[1]
+     r2    <- s$r.squared
+     
+     results <- rbind(results, data.frame(
+         channel       = ch,
+         slope         = slope,
+         slope_SE      = s$coefficients[2, 2],
+         intercept     = int,
+         intercept_SE  = s$coefficients[1, 2],
+         R_squared     = r2,
+         regression_SE = s$sigma,
+         row.names     = NULL))
+     
+     # Equation label, e.g.  y = -4.882x + 65.82,  R^2 = 0.936
+     eq <- sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+                   slope, ifelse(int < 0, "-", "+"), abs(int), r2)
+     
+     # Put the equation in an empty corner (bottom-left if the line goes
+     # down, top-left if it goes up)
+     y_pad <- 0.08 * diff(range(dat$y))
+     if (slope < 0) {
+         eq_y <- min(dat$y) - y_pad; eq_vjust <- 0
+     } else {
+         eq_y <- max(dat$y) + y_pad; eq_vjust <- 1
+     }
+     
+     p <- ggplot(dat, aes(x = conc, y = y)) +
+         geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
+                     colour = line_color, linewidth = 0.6) +
+         geom_point(shape = 21, size = 2.8, fill = col,
+                    colour = "white", stroke = 0.8) +
+         annotate("text", x = min(dat$conc), y = eq_y, hjust = 0,
+                  vjust = eq_vjust, size = 2.6, parse = TRUE, label = eq) +
+         scale_x_continuous(breaks = standards$conc_mg_mL) +
+         labs(title = paste0("Calibration curve: ", tolower(ch),
+                             " vs. protein concentration"),
+              x = "Protein standard concentration (mg/mL)",
+              y = y_labels[ch]) +
+         theme_classic(base_size = 8) +
+         theme(plot.title = element_text(size = 8, face = "bold"),
+               axis.text  = element_text(colour = "black"))
+     
+     print(p)
+     file_name <- sprintf("Figure_%d_%s_calibration.png", i, ch)
+     ggsave(file_name, p, width = 3.3, height = 2.8, units = "in", dpi = 600)
+     cat("Saved", file_name, "\n")
+ }
Saved Figure_1_Red_calibration.png 
Saved Figure_2_Green_calibration.png 
Saved Figure_3_Blue_calibration.png 
Saved Figure_4_Hue_calibration.png 
Saved Figure_5_Saturation_calibration.png 
Saved Figure_6_Lightness_calibration.png 
> 
> # ---- 6. SUMMARY TABLE: WHICH CHANNEL IS MOST LINEAR? -------------------
> results <- results[order(-results$R_squared), ]
> cat("\n===== Linear fits for all 6 channels (best R^2 first) =====\n")

===== Linear fits for all 6 channels (best R^2 first) =====
> print(results, row.names = FALSE, digits = 4)
    channel   slope slope_SE intercept intercept_SE R_squared regression_SE
  Lightness  -4.882   0.7371     65.82        2.445    0.9360         2.331
      Green -11.100   2.2591    194.70        7.492    0.8895         7.144
        Red -17.400   3.8884    258.60       12.897    0.8697        12.296
       Blue  -7.500   2.0551     77.10        6.816    0.8162         6.499
 Saturation  -4.618   2.7456     83.05        9.106    0.4853         8.682
        Hue   1.467   1.4661     38.07        4.862    0.2503         4.636
> write.csv(results, "calibration_fits_all_channels.csv", row.names = FALSE)
> 
> cat("\nBest channel (highest R^2):", results$channel[1],
+     sprintf("(R^2 = %.3f)", results$R_squared[1]), "\n")

Best channel (highest R^2): Lightness (R^2 = 0.936) 
> cat("\nAll files saved to:", getwd(), "\n")

All files saved to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # Milk samples vs. the LIGHTNESS calibration curve
> #
> # What this script does:
> #   1. Builds the lightness calibration curve from the 5 protein standards
> #   2. Converts each milk disk's color to lightness
> #   3. Uses the curve to calculate protein (mg/mL) in each milk sample
> #   4. Compares each milk to the 3 mg/mL it was diluted to (% recovery)
> #   5. Runs t-tests (each milk vs. 3 mg/mL, and whole vs. almond)
> #   6. Makes 2 graphs:
> #        Figure_milk_on_calibration_curve.png
> #        Figure_milk_comparison_bars.png
> #
> # HOW TO USE: open in R, click "Source". Files are saved to the folder
> # printed at the end.
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> 
> # Protein standards (80:20 casein:BSA), one Color Grab reading each
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # Milk samples: average color of the center of each indicator disk
> # (pixels averaged over a circle inside each disk, because the disks
> # were blotchy and a single tap was not representative)
> milk <- data.frame(
+     sample    = c("Whole milk", "Whole milk", "Whole milk",
+                   "Almond milk", "Almond milk", "Almond milk"),
+     replicate = c(1, 2, 3, 1, 2, 3),
+     hex       = c("#697159", "#5E6B5C", "#7A7D74",    # whole milk disks
+                   "#AF8B53", "#536456", "#A09A88")    # almond milk disks
+ )
> 
> # Both milks were diluted to this protein concentration (based on label)
> target_mg_mL <- 3
> 
> # OPTIONAL: protein on each carton label, to calculate the protein in the
> # original (undiluted) milk. Leave as NA to skip.
> # Example for 8 g per 1 cup:  label_g = 8, serving_mL = 240
> label <- data.frame(
+     sample     = c("Whole milk", "Almond milk"),
+     label_g    = c(NA, NA),
+     serving_mL = c(NA, NA)
+ )
> 
> # Graph colors
> milk_colors <- c("Whole milk" = "#2A78D6",   # blue
+                  "Almond milk" = "#EB6834")  # orange
> 
> # ---- 2. HEX CODE -> LIGHTNESS (%) --------------------------------------
> # Lightness = average of the brightest and darkest of R, G, B
> lightness <- function(hex) {
+     rgb_vals <- col2rgb(hex) / 255
+     (apply(rgb_vals, 2, max) + apply(rgb_vals, 2, min)) / 2 * 100
+ }
> standards$L <- lightness(standards$hex)
> milk$L      <- lightness(milk$hex)
> 
> # ---- 3. CALIBRATION CURVE (lightness vs. protein) -----------------------
> fit       <- lm(L ~ conc_mg_mL, data = standards)
> slope     <- coef(fit)[2]
> intercept <- coef(fit)[1]
> r2        <- summary(fit)$r.squared
> cat("\n===== Lightness calibration curve =====\n")

===== Lightness calibration curve =====
> cat(sprintf("Lightness = %.3f * concentration + %.2f   (R^2 = %.3f)\n",
+             slope, intercept, r2))
Lightness = -4.882 * concentration + 65.82   (R^2 = 0.936)
> 
> # ---- 4. PROTEIN IN EACH MILK SAMPLE -------------------------------------
> # Rearranged line:  concentration = (lightness - intercept) / slope
> milk$conc_mg_mL <- (milk$L - intercept) / slope
> milk$recovery_percent <- milk$conc_mg_mL / target_mg_mL * 100
> milk$in_standard_range <- milk$conc_mg_mL >= min(standards$conc_mg_mL) &
+     milk$conc_mg_mL <= max(standards$conc_mg_mL)
> 
> cat("\n===== Each milk disk =====\n")

===== Each milk disk =====
> print(milk, row.names = FALSE, digits = 4)
      sample replicate     hex     L conc_mg_mL recovery_percent in_standard_range
  Whole milk         1 #697159 39.61      5.369           178.98             FALSE
  Whole milk         2 #5E6B5C 39.02      5.490           183.00             FALSE
  Whole milk         3 #7A7D74 47.25      3.803           126.77              TRUE
 Almond milk         1 #AF8B53 50.59      3.120           104.02              TRUE
 Almond milk         2 #536456 35.88      6.133           204.42             FALSE
 Almond milk         3 #A09A88 58.04      1.594            53.15              TRUE
> 
> # Summary for each milk: mean, SD, 95% confidence interval, % recovery
> summ <- do.call(rbind, lapply(split(milk, milk$sample), function(d) {
+     n  <- nrow(d); m <- mean(d$conc_mg_mL); s <- sd(d$conc_mg_mL)
+     half <- qt(0.975, n - 1) * s / sqrt(n)
+     tt <- t.test(d$conc_mg_mL, mu = target_mg_mL)   # vs. 3 mg/mL
+     data.frame(sample = d$sample[1], n = n,
+                mean_mg_mL = m, sd_mg_mL = s,
+                CI95_low = m - half, CI95_high = m + half,
+                mean_L = mean(d$L),
+                recovery_percent = m / target_mg_mL * 100,
+                t_vs_target = unname(tt$statistic),
+                p_vs_target = tt$p.value)
+ }))
> summ$sample <- factor(summ$sample, levels = names(milk_colors))
> summ <- summ[order(summ$sample), ]
> 
> cat("\n===== Summary (protein in the DILUTED samples) =====\n")

===== Summary (protein in the DILUTED samples) =====
> print(summ, row.names = FALSE, digits = 4)
      sample n mean_mg_mL sd_mg_mL CI95_low CI95_high mean_L recovery_percent t_vs_target p_vs_target
  Whole milk 3      4.888    0.941    2.550     7.225  41.96            162.9      3.4743     0.07379
 Almond milk 3      3.616    2.309   -2.121     9.352  48.17            120.5      0.4619     0.68954
> 
> # Whole vs. almond (Welch two-sample t-test)
> tt2 <- t.test(conc_mg_mL ~ sample, data = milk)
> cat("\n===== Whole milk vs. almond milk (Welch t-test) =====\n")

===== Whole milk vs. almond milk (Welch t-test) =====
> print(tt2)

	Welch Two Sample t-test

data:  conc_mg_mL by sample
t = -0.88335, df = 2.6464, p-value = 0.4499
alternative hypothesis: true difference in means between group Almond milk and group Whole milk is not equal to 0
95 percent confidence interval:
 -6.220138  3.676631
sample estimates:
mean in group Almond milk  mean in group Whole milk 
                 3.615797                  4.887550 

> 
> # OPTIONAL: protein in the original, undiluted milk
> summ <- merge(summ, label, by = "sample", sort = FALSE)
> summ$label_mg_mL     <- summ$label_g * 1000 / summ$serving_mL
> summ$dilution_factor <- summ$label_mg_mL / target_mg_mL
> summ$original_mg_mL  <- summ$mean_mg_mL * summ$dilution_factor
> if (any(!is.na(summ$original_mg_mL))) {
+     cat("\n===== Protein in the ORIGINAL (undiluted) milk =====\n")
+     print(summ[, c("sample", "label_mg_mL", "dilution_factor",
+                    "mean_mg_mL", "original_mg_mL")], row.names = FALSE, digits = 4)
+ }
> 
> write.csv(milk, "milk_results_each_disk.csv", row.names = FALSE)
> write.csv(summ, "milk_results_summary.csv",   row.names = FALSE)
> 
> # ---- 5. GRAPH 1: MILK SAMPLES ON THE CALIBRATION CURVE ------------------
> std_plot  <- data.frame(conc = standards$conc_mg_mL, L = standards$L,
+                         group = "Protein standards")
> milk_plot <- data.frame(conc = summ$mean_mg_mL, L = summ$mean_L,
+                         sd = summ$sd_mg_mL, group = as.character(summ$sample))
> all_colors <- c("Protein standards" = "black", milk_colors)
> 
> eq <- sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+               slope, ifelse(intercept < 0, "-", "+"), abs(intercept), r2)
> 
> fig_curve <- ggplot() +
+     # target concentration
+     annotate("segment", x = target_mg_mL, xend = target_mg_mL,
+              y = -Inf, yend = max(std_plot$L) - 1, linetype = "dashed",
+              colour = "grey55", linewidth = 0.4) +
+     annotate("text", x = target_mg_mL, y = min(std_plot$L) - 4,
+              label = "target 3 mg/mL", hjust = -0.05, vjust = 0, size = 2.3,
+              colour = "grey35") +
+     # calibration line
+     geom_smooth(data = std_plot, aes(x = conc, y = L), method = "lm",
+                 formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
+     # milk error bars (SD of calculated concentration)
+     geom_errorbar(data = milk_plot,
+                   aes(y = L, xmin = conc - sd, xmax = conc + sd, colour = group),
+                   orientation = "y", width = 0.8, linewidth = 0.5,
+                   show.legend = FALSE) +
+     # standards and milk points
+     geom_point(data = std_plot, aes(x = conc, y = L, fill = group),
+                shape = 21, size = 2.6, colour = "white", stroke = 0.8) +
+     geom_point(data = milk_plot, aes(x = conc, y = L, fill = group),
+                shape = 21, size = 3.2, colour = "white", stroke = 0.8) +
+     scale_fill_manual(values = all_colors, breaks = names(all_colors), name = NULL) +
+     scale_colour_manual(values = all_colors) +
+     annotate("text", x = 6.3, y = max(std_plot$L) + 1.5, hjust = 1, vjust = 1,
+              size = 2.5, parse = TRUE, label = eq) +
+     scale_x_continuous(breaks = 1:6, limits = c(0.3, 6.3)) +
+     labs(title = "Milk samples on the lightness calibration curve",
+          x = "Protein concentration (mg/mL)",
+          y = "Lightness (%)") +
+     theme_classic(base_size = 8) +
+     theme(plot.title = element_text(size = 8, face = "bold"),
+           axis.text = element_text(colour = "black"),
+           legend.position = "bottom",
+           legend.margin = margin(0, 0, 0, 0))
> 
> print(fig_curve)
> ggsave("Figure_milk_on_calibration_curve.png", fig_curve,
+        width = 3.3, height = 3.1, units = "in", dpi = 600)
> 
> # ---- 6. GRAPH 2: WHOLE vs. ALMOND BAR CHART -----------------------------
> milk$sample <- factor(milk$sample, levels = names(milk_colors))
> p_label <- sprintf(paste0("Dashed line: 3 mg/mL target. Error bars: SD (n = 3).\n",
+                           "Whole vs. almond: p = %.2f (Welch t-test)"), tt2$p.value)
> 
> fig_bars <- ggplot(summ, aes(x = sample, y = mean_mg_mL)) +
+     geom_col(aes(fill = sample), width = 0.55) +
+     geom_errorbar(aes(ymin = mean_mg_mL - sd_mg_mL, ymax = mean_mg_mL + sd_mg_mL),
+                   width = 0.15, linewidth = 0.5) +
+     geom_point(data = milk, aes(x = sample, y = conc_mg_mL),
+                position = position_jitter(width = 0.08, height = 0, seed = 1),
+                shape = 21, size = 1.8, fill = "white", colour = "black",
+                stroke = 0.5) +
+     geom_hline(yintercept = target_mg_mL, linetype = "dashed",
+                colour = "grey35", linewidth = 0.4) +
+     scale_fill_manual(values = milk_colors, guide = "none") +
+     scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
+     labs(title = "Protein measured in diluted milk samples",
+          subtitle = p_label,
+          x = NULL,
+          y = "Protein in diluted sample (mg/mL)") +
+     theme_classic(base_size = 8) +
+     theme(plot.title = element_text(size = 8, face = "bold"),
+           plot.subtitle = element_text(size = 7, colour = "grey25"),
+           axis.text = element_text(colour = "black"))
> 
> print(fig_bars)
> ggsave("Figure_milk_comparison_bars.png", fig_bars,
+        width = 3.3, height = 3.0, units = "in", dpi = 600)
> 
> cat("\nAll files saved to:", getwd(), "\n")

All files saved to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # Milk samples vs. the LIGHTNESS calibration curve
> #
> # What this script does:
> #   1. Builds the lightness calibration curve from the 5 protein standards
> #   2. Converts each milk disk's color to lightness
> #   3. Uses the curve to calculate protein (mg/mL) in each milk sample
> #   4. Compares each milk to the 3 mg/mL it was diluted to (% recovery)
> #   5. Runs t-tests (each milk vs. 3 mg/mL, and whole vs. almond)
> #   6. Makes 2 graphs:
> #        Figure_milk_on_calibration_curve.png
> #        Figure_milk_comparison_bars.png
> #
> # HOW TO USE: open in R, click "Source". Files are saved to the folder
> # printed at the end.
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> 
> # Protein standards (80:20 casein:BSA), one Color Grab reading each
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # Milk samples: average color of the center of each indicator disk
> # (pixels averaged over a circle inside each disk, because the disks
> # were blotchy and a single tap was not representative)
> milk <- data.frame(
+     sample    = c("Whole milk", "Whole milk", "Whole milk",
+                   "Almond milk", "Almond milk", "Almond milk"),
+     replicate = c(1, 2, 3, 1, 2, 3),
+     hex       = c("#697159", "#5E6B5C", "#7A7D74",    # whole milk disks
+                   "#AF8B53", "#536456", "#A09A88")    # almond milk disks
+ )
> 
> # Both milks were diluted to this protein concentration (based on label)
> target_mg_mL <- 3
> 
> # OPTIONAL: protein on each carton label, to calculate the protein in the
> # original (undiluted) milk. Leave as NA to skip.
> # Example for 8 g per 1 cup:  label_g = 8, serving_mL = 240
> label <- data.frame(
+     sample     = c("Whole milk", "Almond milk"),
+     label_g    = c(NA, NA),
+     serving_mL = c(NA, NA)
+ )
> 
> # Graph colors
> milk_colors <- c("Whole milk" = "#2A78D6",   # blue
+                  "Almond milk" = "#EB6834")  # orange
> 
> # ---- 2. HEX CODE -> LIGHTNESS (%) --------------------------------------
> # Lightness = average of the brightest and darkest of R, G, B
> lightness <- function(hex) {
+     rgb_vals <- col2rgb(hex) / 255
+     (apply(rgb_vals, 2, max) + apply(rgb_vals, 2, min)) / 2 * 100
+ }
> standards$L <- lightness(standards$hex)
> milk$L      <- lightness(milk$hex)
> 
> # ---- 3. CALIBRATION CURVE (lightness vs. protein) -----------------------
> fit       <- lm(L ~ conc_mg_mL, data = standards)
> slope     <- coef(fit)[2]
> intercept <- coef(fit)[1]
> r2        <- summary(fit)$r.squared
> cat("\n===== Lightness calibration curve =====\n")

===== Lightness calibration curve =====
> cat(sprintf("Lightness = %.3f * concentration + %.2f   (R^2 = %.3f)\n",
+             slope, intercept, r2))
Lightness = -4.882 * concentration + 65.82   (R^2 = 0.936)
> 
> # ---- 4. PROTEIN IN EACH MILK SAMPLE -------------------------------------
> # Rearranged line:  concentration = (lightness - intercept) / slope
> milk$conc_mg_mL <- (milk$L - intercept) / slope
> milk$recovery_percent <- milk$conc_mg_mL / target_mg_mL * 100
> milk$in_standard_range <- milk$conc_mg_mL >= min(standards$conc_mg_mL) &
+     milk$conc_mg_mL <= max(standards$conc_mg_mL)
> 
> cat("\n===== Each milk disk =====\n")

===== Each milk disk =====
> print(milk, row.names = FALSE, digits = 4)
      sample replicate     hex     L conc_mg_mL recovery_percent in_standard_range
  Whole milk         1 #697159 39.61      5.369           178.98             FALSE
  Whole milk         2 #5E6B5C 39.02      5.490           183.00             FALSE
  Whole milk         3 #7A7D74 47.25      3.803           126.77              TRUE
 Almond milk         1 #AF8B53 50.59      3.120           104.02              TRUE
 Almond milk         2 #536456 35.88      6.133           204.42             FALSE
 Almond milk         3 #A09A88 58.04      1.594            53.15              TRUE
> 
> # Summary for each milk: mean, SD, 95% confidence interval, % recovery
> summ <- do.call(rbind, lapply(split(milk, milk$sample), function(d) {
+     n  <- nrow(d); m <- mean(d$conc_mg_mL); s <- sd(d$conc_mg_mL)
+     half <- qt(0.975, n - 1) * s / sqrt(n)
+     tt <- t.test(d$conc_mg_mL, mu = target_mg_mL)   # vs. 3 mg/mL
+     data.frame(sample = d$sample[1], n = n,
+                mean_mg_mL = m, sd_mg_mL = s,
+                CI95_low = m - half, CI95_high = m + half,
+                mean_L = mean(d$L),
+                recovery_percent = m / target_mg_mL * 100,
+                t_vs_target = unname(tt$statistic),
+                p_vs_target = tt$p.value)
+ }))
> summ$sample <- factor(summ$sample, levels = names(milk_colors))
> summ <- summ[order(summ$sample), ]
> 
> cat("\n===== Summary (protein in the DILUTED samples) =====\n")

===== Summary (protein in the DILUTED samples) =====
> print(summ, row.names = FALSE, digits = 4)
      sample n mean_mg_mL sd_mg_mL CI95_low CI95_high mean_L recovery_percent t_vs_target p_vs_target
  Whole milk 3      4.888    0.941    2.550     7.225  41.96            162.9      3.4743     0.07379
 Almond milk 3      3.616    2.309   -2.121     9.352  48.17            120.5      0.4619     0.68954
> 
> # Whole vs. almond (Welch two-sample t-test)
> tt2 <- t.test(conc_mg_mL ~ sample, data = milk)
> cat("\n===== Whole milk vs. almond milk (Welch t-test) =====\n")

===== Whole milk vs. almond milk (Welch t-test) =====
> print(tt2)

	Welch Two Sample t-test

data:  conc_mg_mL by sample
t = -0.88335, df = 2.6464, p-value = 0.4499
alternative hypothesis: true difference in means between group Almond milk and group Whole milk is not equal to 0
95 percent confidence interval:
 -6.220138  3.676631
sample estimates:
mean in group Almond milk  mean in group Whole milk 
                 3.615797                  4.887550 

> 
> # OPTIONAL: protein in the original, undiluted milk
> summ <- merge(summ, label, by = "sample", sort = FALSE)
> summ$label_mg_mL     <- summ$label_g * 1000 / summ$serving_mL
> summ$dilution_factor <- summ$label_mg_mL / target_mg_mL
> summ$original_mg_mL  <- summ$mean_mg_mL * summ$dilution_factor
> if (any(!is.na(summ$original_mg_mL))) {
+     cat("\n===== Protein in the ORIGINAL (undiluted) milk =====\n")
+     print(summ[, c("sample", "label_mg_mL", "dilution_factor",
+                    "mean_mg_mL", "original_mg_mL")], row.names = FALSE, digits = 4)
+ }
> 
> write.csv(milk, "milk_results_each_disk.csv", row.names = FALSE)
> write.csv(summ, "milk_results_summary.csv",   row.names = FALSE)
> 
> # ---- 5. GRAPH 1: MILK SAMPLES ON THE CALIBRATION CURVE ------------------
> std_plot  <- data.frame(conc = standards$conc_mg_mL, L = standards$L,
+                         group = "Protein standards")
> milk_plot <- data.frame(conc = summ$mean_mg_mL, L = summ$mean_L,
+                         sd = summ$sd_mg_mL, group = as.character(summ$sample))
> all_colors <- c("Protein standards" = "black", milk_colors)
> 
> eq <- sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+               slope, ifelse(intercept < 0, "-", "+"), abs(intercept), r2)
> 
> fig_curve <- ggplot() +
+     # target concentration
+     annotate("segment", x = target_mg_mL, xend = target_mg_mL,
+              y = -Inf, yend = max(std_plot$L) - 1, linetype = "dashed",
+              colour = "grey55", linewidth = 0.4) +
+     annotate("text", x = target_mg_mL, y = min(std_plot$L) - 4,
+              label = "target 3 mg/mL", hjust = -0.05, vjust = 0, size = 2.3,
+              colour = "grey35") +
+     # calibration line
+     geom_smooth(data = std_plot, aes(x = conc, y = L), method = "lm",
+                 formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
+     # milk error bars (SD of calculated concentration)
+     geom_errorbar(data = milk_plot,
+                   aes(y = L, xmin = conc - sd, xmax = conc + sd, colour = group),
+                   orientation = "y", width = 0.8, linewidth = 0.5,
+                   show.legend = FALSE) +
+     # standards and milk points
+     geom_point(data = std_plot, aes(x = conc, y = L, fill = group),
+                shape = 21, size = 2.6, colour = "white", stroke = 0.8) +
+     geom_point(data = milk_plot, aes(x = conc, y = L, fill = group),
+                shape = 21, size = 3.2, colour = "white", stroke = 0.8) +
+     scale_fill_manual(values = all_colors, breaks = names(all_colors), name = NULL) +
+     scale_colour_manual(values = all_colors) +
+     annotate("text", x = 6.3, y = max(std_plot$L) + 1.5, hjust = 1, vjust = 1,
+              size = 2.5, parse = TRUE, label = eq) +
+     scale_x_continuous(breaks = 1:6, limits = c(0.3, 6.3)) +
+     labs(title = "Milk samples on the lightness calibration curve",
+          x = "Protein concentration (mg/mL)",
+          y = "Lightness (%)") +
+     theme_classic(base_size = 8) +
+     theme(plot.title = element_text(size = 8, face = "bold"),
+           axis.text = element_text(colour = "black"),
+           legend.position = "bottom",
+           legend.margin = margin(0, 0, 0, 0))
> 
> print(fig_curve)
> ggsave("Figure_milk_on_calibration_curve.png", fig_curve,
+        width = 3.3, height = 3.1, units = "in", dpi = 600)
> 
> # ---- 6. GRAPH 2: WHOLE vs. ALMOND BAR CHART -----------------------------
> milk$sample <- factor(milk$sample, levels = names(milk_colors))
> 
> # Significance label for whole vs. almond (Welch t-test)
> p_val <- tt2$p.value
> sig_label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else
+     if (p_val < 0.05) "*" else "n.s."
> sig_text  <- sprintf("%s  (p = %s)", sig_label,
+                      ifelse(p_val < 0.001, "< 0.001", sprintf("%.2f", p_val)))
> cat("\nWhole vs. almond:", sig_text, "\n")

Whole vs. almond: n.s.  (p = 0.45) 
> 
> # Where to draw the bracket: just above the tallest bar, error bar, or dot
> top       <- max(summ$mean_mg_mL + summ$sd_mg_mL, milk$conc_mg_mL)
> bracket_y <- top * 1.08
> tick      <- top * 0.03
> 
> p_label <- paste0("Dashed line: 3 mg/mL target. Error bars: SD (n = 3).\n",
+                   "Welch t-test: n.s. = not significant (p > 0.05)\n",
+                   "* p < 0.05, ** p < 0.01, *** p < 0.001")
> 
> fig_bars <- ggplot(summ, aes(x = sample, y = mean_mg_mL)) +
+     geom_col(aes(fill = sample), width = 0.55) +
+     geom_errorbar(aes(ymin = mean_mg_mL - sd_mg_mL, ymax = mean_mg_mL + sd_mg_mL),
+                   width = 0.15, linewidth = 0.5) +
+     geom_point(data = milk, aes(x = sample, y = conc_mg_mL),
+                position = position_jitter(width = 0.08, height = 0, seed = 1),
+                shape = 21, size = 1.8, fill = "white", colour = "black",
+                stroke = 0.5) +
+     geom_hline(yintercept = target_mg_mL, linetype = "dashed",
+                colour = "grey35", linewidth = 0.4) +
+     # significance bracket between the two bars
+     annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
+              linewidth = 0.4) +
+     annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - tick,
+              linewidth = 0.4) +
+     annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - tick,
+              linewidth = 0.4) +
+     annotate("text", x = 1.5, y = bracket_y, label = sig_text,
+              vjust = -0.5, size = 2.6) +
+     scale_fill_manual(values = milk_colors, guide = "none") +
+     scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
+                        limits = c(0, bracket_y * 1.12)) +
+     labs(title = "Protein measured in diluted milk samples",
+          subtitle = p_label,
+          x = NULL,
+          y = "Protein in diluted sample (mg/mL)") +
+     theme_classic(base_size = 8) +
+     theme(plot.title = element_text(size = 8, face = "bold"),
+           plot.subtitle = element_text(size = 7, colour = "grey25"),
+           axis.text = element_text(colour = "black"))
> 
> print(fig_bars)
> ggsave("Figure_milk_comparison_bars.png", fig_bars,
+        width = 3.3, height = 3.3, units = "in", dpi = 600)
> 
> cat("\nAll files saved to:", getwd(), "\n")

All files saved to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # ONE combined figure with all 6 calibration curves (panels A-F):
> #   A Red, B Green, C Blue, D Hue, E Saturation, F Lightness
> #
> # Style: black line of best fit on every panel; red, green, and blue
> # dots on the R, G, B panels; black dots on the H, S, L panels.
> # Each panel shows its equation and R^2.
> #
> # Data: ONE Color Grab reading per protein standard (original readings).
> # Standards: 1.0, 2.0, 3.0, 4.0 mg/mL dilutions + 5.0 mg/mL stock
> #            (80:20 casein:BSA, from the lab procedure).
> #
> # Output (saved to the folder printed at the end):
> #   Figure_six_channel_calibration.png   (ACS double-column width, 7 in)
> #   calibration_fits_all_channels.csv     (numbers for Table 1)
> #
> # HOW TO USE: open in R, click "Source".
> # =====================================================================
> 
> library(ggplot2)
> 
> # ---- 1. DATA ----------------------------------------------------------
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # ---- 2. COLORS (change any of these if you want) ------------------------
> dot_colors <- c(Red        = "#E34948",   # red
+                 Green      = "#008300",   # green
+                 Blue       = "#2A78D6",   # blue
+                 Hue        = "black",
+                 Saturation = "black",
+                 Lightness  = "black")
> line_color <- "black"
> 
> # ---- 3. HEX CODE -> R, G, B, HUE, SATURATION, LIGHTNESS ---------------
> # R, G, B: 0-255.  Hue: degrees (0-360).  Saturation, Lightness: percent.
> rgb_vals <- t(col2rgb(standards$hex))
> r <- rgb_vals[, 1] / 255
> g <- rgb_vals[, 2] / 255
> b <- rgb_vals[, 3] / 255
> cmax  <- pmax(r, g, b)
> cmin  <- pmin(r, g, b)
> delta <- cmax - cmin
> 
> L <- (cmax + cmin) / 2
> S <- ifelse(delta == 0, 0, delta / (1 - abs(2 * L - 1)))
> H <- ifelse(delta == 0, 0,
+             ifelse(cmax == r, 60 * (((g - b) / delta) %% 6),
+                    ifelse(cmax == g, 60 * ((b - r) / delta + 2),
+                           60 * ((r - g) / delta + 4))))
> 
> values <- list(Red        = rgb_vals[, 1],
+                Green      = rgb_vals[, 2],
+                Blue       = rgb_vals[, 3],
+                Hue        = H,
+                Saturation = S * 100,
+                Lightness  = L * 100)
> 
> # ---- 4. PANEL NAMES (letter + what the y-axis means) --------------------
> channels <- c("Red", "Green", "Blue", "Hue", "Saturation", "Lightness")
> panel_names <- c(Red        = "A. Red intensity (0-255)",
+                  Green      = "B. Green intensity (0-255)",
+                  Blue       = "C. Blue intensity (0-255)",
+                  Hue        = "D. Hue (degrees)",
+                  Saturation = "E. Saturation (%)",
+                  Lightness  = "F. Lightness (%)")
> 
> # One long table: every standard x every channel
> long <- do.call(rbind, lapply(channels, function(ch)
+     data.frame(channel = ch, conc = standards$conc_mg_mL, y = values[[ch]])))
> long$panel <- factor(panel_names[long$channel], levels = panel_names)
> 
> # ---- 5. LINEAR FIT FOR EACH CHANNEL -----------------------------------
> fits <- do.call(rbind, lapply(channels, function(ch) {
+     d   <- subset(long, channel == ch)
+     fit <- lm(y ~ conc, data = d)
+     s   <- summary(fit)
+     data.frame(channel       = ch,
+                slope         = coef(fit)[2],
+                slope_SE      = s$coefficients[2, 2],
+                intercept     = coef(fit)[1],
+                intercept_SE  = s$coefficients[1, 2],
+                R_squared     = s$r.squared,
+                regression_SE = s$sigma,
+                row.names     = NULL)
+ }))
> 
> cat("\n===== Linear fits for all 6 channels (Table 1) =====\n")

===== Linear fits for all 6 channels (Table 1) =====
> print(fits, row.names = FALSE, digits = 4)
    channel   slope slope_SE intercept intercept_SE R_squared regression_SE
        Red -17.400   3.8884    258.60       12.897    0.8697        12.296
      Green -11.100   2.2591    194.70        7.492    0.8895         7.144
       Blue  -7.500   2.0551     77.10        6.816    0.8162         6.499
        Hue   1.467   1.4661     38.07        4.862    0.2503         4.636
 Saturation  -4.618   2.7456     83.05        9.106    0.4853         8.682
  Lightness  -4.882   0.7371     65.82        2.445    0.9360         2.331
> write.csv(fits, "calibration_fits_all_channels.csv", row.names = FALSE)
> 
> # Equation + R^2 label for each panel, placed in an empty corner
> # (bottom-left if the line goes down, top-left if it goes up)
> labels <- do.call(rbind, lapply(channels, function(ch) {
+     f   <- fits[fits$channel == ch, ]
+     d   <- subset(long, channel == ch)
+     pad <- 0.15 * diff(range(d$y))
+     down <- f$slope < 0
+     data.frame(
+         panel = factor(panel_names[ch], levels = panel_names),
+         x     = min(d$conc),
+         y     = if (down) min(d$y) - pad else max(d$y) + pad,
+         vjust = if (down) 0 else 1,
+         label = sprintf('"y = %.3fx %s %.2f,"~~italic(R)^2~"= %.3f"',
+                         f$slope, ifelse(f$intercept < 0, "-", "+"),
+                         abs(f$intercept), f$R_squared))
+ }))
> 
> # ---- 6. THE COMBINED FIGURE --------------------------------------------
> fig <- ggplot(long, aes(x = conc, y = y)) +
+     geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
+                 colour = line_color, linewidth = 0.5) +
+     geom_point(aes(fill = channel), shape = 21, size = 2.2,
+                colour = "white", stroke = 0.6) +
+     geom_text(data = labels,
+               aes(x = x, y = y, label = label, vjust = vjust),
+               hjust = 0, size = 2.3, parse = TRUE) +
+     scale_fill_manual(values = dot_colors, guide = "none") +
+     facet_wrap(~ panel, ncol = 3, scales = "free") +   # axes on every panel
+     scale_x_continuous(breaks = standards$conc_mg_mL) +
+     labs(title = "Calibration curves for each color channel vs. protein concentration",
+          x = "Protein standard concentration (mg/mL)",
+          y = "Channel value") +
+     theme_classic(base_size = 8) +
+     theme(plot.title       = element_text(size = 9, face = "bold"),
+           strip.background = element_blank(),
+           strip.text       = element_text(face = "bold", hjust = 0, size = 8),
+           axis.text        = element_text(colour = "black"),
+           panel.spacing.x  = unit(0.6, "lines"),
+           panel.spacing.y  = unit(0.9, "lines"))
> 
> print(fig)
> ggsave("Figure_six_channel_calibration.png", fig,
+        width = 7, height = 4.6, units = "in", dpi = 600)
> 
> cat("\nSaved Figure_six_channel_calibration.png to:", getwd(), "\n")

Saved Figure_six_channel_calibration.png to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> # =====================================================================
> # CHEM 124 Lab 2 - Protein in milk (TBPB paper microfluidic device)
> # ONE combined figure with all 6 calibration curves (panels A-F):
> #   A Red, B Green, C Blue, D Hue, E Saturation, F Lightness
> #
> # Style: black line of best fit on every panel; red, green, and blue
> # dots on the R, G, B panels; black dots on the H, S, L panels.
> # Each panel shows its equation and R^2 on two short lines so nothing
> # gets cut off.
> #
> # Saves TWO versions (use the one that matches where you put it):
> #   Figure_six_channel_1column.png  -> fits ONE column  (3.3 in wide)
> #   Figure_six_channel_2column.png  -> fits the FULL page width (7 in)
> # plus calibration_fits_all_channels.csv (numbers for Table 1)
> #
> # IMPORTANT: insert the saved PNG file into your document
> # (Insert > Image > Upload from computer). Don't copy the graph from the
> # RStudio Plots pane - that window squishes the panels and cuts labels off.
> #
> # HOW TO USE: open in R, click "Source".
> # =====================================================================
> 
> library(ggplot2)
> 
> # Show a title inside the image? Set to FALSE if your document already
> # has the figure title above the image.
> show_title <- TRUE
> 
> # ---- 1. DATA ----------------------------------------------------------
> standards <- data.frame(
+     conc_mg_mL = c(1.0, 2.0, 3.0, 4.0, 5.0),
+     hex        = c("#EABB44",   # 1.0 mg/mL
+                    "#E2AD3D",   # 2.0 mg/mL
+                    "#DF9E37",   # 3.0 mg/mL
+                    "#B28E38",   # 4.0 mg/mL
+                    "#AB9321")   # 5.0 mg/mL
+ )
> 
> # ---- 2. COLORS (change any of these if you want) ------------------------
> dot_colors <- c(Red        = "#E34948",   # red
+                 Green      = "#008300",   # green
+                 Blue       = "#2A78D6",   # blue
+                 Hue        = "black",
+                 Saturation = "black",
+                 Lightness  = "black")
> line_color <- "black"
> 
> # ---- 3. HEX CODE -> R, G, B, HUE, SATURATION, LIGHTNESS ---------------
> # R, G, B: 0-255.  Hue: degrees (0-360).  Saturation, Lightness: percent.
> rgb_vals <- t(col2rgb(standards$hex))
> r <- rgb_vals[, 1] / 255
> g <- rgb_vals[, 2] / 255
> b <- rgb_vals[, 3] / 255
> cmax  <- pmax(r, g, b)
> cmin  <- pmin(r, g, b)
> delta <- cmax - cmin
> 
> L <- (cmax + cmin) / 2
> S <- ifelse(delta == 0, 0, delta / (1 - abs(2 * L - 1)))
> H <- ifelse(delta == 0, 0,
+             ifelse(cmax == r, 60 * (((g - b) / delta) %% 6),
+                    ifelse(cmax == g, 60 * ((b - r) / delta + 2),
+                           60 * ((r - g) / delta + 4))))
> 
> values <- list(Red        = rgb_vals[, 1],
+                Green      = rgb_vals[, 2],
+                Blue       = rgb_vals[, 3],
+                Hue        = H,
+                Saturation = S * 100,
+                Lightness  = L * 100)
> 
> # ---- 4. PANEL NAMES (short, so they fit) --------------------------------
> channels <- c("Red", "Green", "Blue", "Hue", "Saturation", "Lightness")
> panel_names <- c(Red        = "A. Red (0-255)",
+                  Green      = "B. Green (0-255)",
+                  Blue       = "C. Blue (0-255)",
+                  Hue        = "D. Hue (degrees)",
+                  Saturation = "E. Saturation (%)",
+                  Lightness  = "F. Lightness (%)")
> 
> long <- do.call(rbind, lapply(channels, function(ch)
+     data.frame(channel = ch, conc = standards$conc_mg_mL, y = values[[ch]])))
> long$panel <- factor(panel_names[long$channel], levels = panel_names)
> 
> # ---- 5. LINEAR FIT FOR EACH CHANNEL -----------------------------------
> fits <- do.call(rbind, lapply(channels, function(ch) {
+     d   <- subset(long, channel == ch)
+     fit <- lm(y ~ conc, data = d)
+     s   <- summary(fit)
+     data.frame(channel       = ch,
+                slope         = coef(fit)[2],
+                slope_SE      = s$coefficients[2, 2],
+                intercept     = coef(fit)[1],
+                intercept_SE  = s$coefficients[1, 2],
+                R_squared     = s$r.squared,
+                regression_SE = s$sigma,
+                row.names     = NULL)
+ }))
> 
> cat("\n===== Linear fits for all 6 channels (Table 1) =====\n")

===== Linear fits for all 6 channels (Table 1) =====
> print(fits, row.names = FALSE, digits = 4)
    channel   slope slope_SE intercept intercept_SE R_squared regression_SE
        Red -17.400   3.8884    258.60       12.897    0.8697        12.296
      Green -11.100   2.2591    194.70        7.492    0.8895         7.144
       Blue  -7.500   2.0551     77.10        6.816    0.8162         6.499
        Hue   1.467   1.4661     38.07        4.862    0.2503         4.636
 Saturation  -4.618   2.7456     83.05        9.106    0.4853         8.682
  Lightness  -4.882   0.7371     65.82        2.445    0.9360         2.331
> write.csv(fits, "calibration_fits_all_channels.csv", row.names = FALSE)
> 
> # Equation (line 1) and R^2 (line 2) for each panel, placed in an empty
> # corner: bottom-left if the line goes down, top-left if it goes up.
> labels <- do.call(rbind, lapply(channels, function(ch) {
+     f    <- fits[fits$channel == ch, ]
+     d    <- subset(long, channel == ch)
+     rng  <- diff(range(d$y))
+     down <- f$slope < 0
+     step <- 0.13 * rng                       # gap between the two lines
+     if (down) {
+         y_r2 <- min(d$y) - 0.30 * rng          # R^2 on the bottom line
+         y_eq <- y_r2 + step                    # equation just above it
+     } else {
+         y_eq <- max(d$y) + 0.30 * rng          # equation on the top line
+         y_r2 <- y_eq - step                    # R^2 just below it
+     }
+     data.frame(
+         panel = factor(panel_names[ch], levels = panel_names),
+         x     = min(d$conc),
+         y_eq  = y_eq,
+         y_r2  = y_r2,
+         eq    = sprintf("y = %.3fx %s %.2f", f$slope,
+                         ifelse(f$intercept < 0, "-", "+"), abs(f$intercept)),
+         r2    = sprintf('italic(R)^2~"= %.3f"', f$R_squared))
+ }))
> 
> # ---- 6. BUILD THE FIGURE ------------------------------------------------
> make_figure <- function(n_col, base, text_size, title) {
+     ggplot(long, aes(x = conc, y = y)) +
+         geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
+                     colour = line_color, linewidth = 0.45) +
+         geom_point(aes(fill = channel), shape = 21, size = base * 0.27,
+                    colour = "white", stroke = 0.5) +
+         geom_text(data = labels, aes(x = x, y = y_eq, label = eq),
+                   hjust = 0, vjust = 0.5, size = text_size) +
+         geom_text(data = labels, aes(x = x, y = y_r2, label = r2),
+                   hjust = 0, vjust = 0.5, size = text_size, parse = TRUE) +
+         scale_fill_manual(values = dot_colors, guide = "none") +
+         facet_wrap(~ panel, ncol = n_col, scales = "free") +
+         scale_x_continuous(breaks = standards$conc_mg_mL,
+                            expand = expansion(mult = c(0.06, 0.06))) +
+         labs(title = if (show_title) title else NULL,
+              x = "Protein standard concentration (mg/mL)",
+              y = "Channel value") +
+         theme_classic(base_size = base) +
+         theme(plot.title       = element_text(size = base + 1, face = "bold"),
+               strip.background = element_blank(),
+               strip.text       = element_text(face = "bold", hjust = 0,
+                                               size = base),
+               axis.text        = element_text(colour = "black",
+                                               size = base - 1),
+               panel.spacing.x  = unit(0.8, "lines"),
+               panel.spacing.y  = unit(0.8, "lines"),
+               plot.margin      = margin(4, 8, 4, 4))
+ }
> 
> # ONE-column version: 2 panels across, 3 rows, 3.3 in wide
> fig_1col <- make_figure(n_col = 2, base = 7, text_size = 2.0,
+                         title = "Calibration curves for each color\nchannel vs. protein concentration")
> ggsave("Figure_six_channel_1column.png", fig_1col,
+        width = 3.3, height = 5.2, units = "in", dpi = 600)
> 
> # TWO-column (full page width) version: 3 panels across, 2 rows, 7 in wide
> fig_2col <- make_figure(n_col = 3, base = 8, text_size = 2.5,
+                         title = "Calibration curves for each color channel vs. protein concentration")
> ggsave("Figure_six_channel_2column.png", fig_2col,
+        width = 7, height = 4.6, units = "in", dpi = 600)
> 
> print(fig_2col)
> cat("\nSaved Figure_six_channel_1column.png and Figure_six_channel_2column.png to:",
+     getwd(), "\n")

Saved Figure_six_channel_1column.png and Figure_six_channel_2column.png to: /Users/johnmicheletos/Lab 2 Calibration Curve 
> 
