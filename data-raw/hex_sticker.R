## Hex sticker for dcvar
## Inspired by https://github.com/GuangchuangYu/hexSticker
## Built with ggplot2 only (no system deps required)

library(ggplot2)

# --- Hex geometry ---
hex_polygon <- function(x_center = 0, y_center = 0, radius = 1) {
  angles <- seq(30, 390, by = 60) * pi / 180
  data.frame(
    x = x_center + radius * cos(angles),
    y = y_center + radius * sin(angles)
  )
}

# For a given y, find the x-extent of the hex (left and right boundary)
hex_x_at_y <- function(y, radius = 0.93) {
  # Hex with flat top/bottom: vertices at 30,90,150,210,270,330 degrees
  # The hex boundary x at a given y is piecewise linear between vertices.
  # For a pointy-left/right hex (our orientation), the width narrows at top/bottom.
  half_w <- radius * sqrt(3) / 2
  half_h <- radius
  # The hex sides run from y = -half_h to y = half_h

  # At y=0 the full width is 2*half_w

  # At y = +/-half_h the width is 0
  # Linear interpolation between the flat section and the points
  ay <- abs(y)
  if (ay > half_h) return(c(0, 0))
  if (ay <= half_h / 2) {
    # In the wide middle band
    xmax <- half_w
  } else {
    # Tapering towards the top/bottom vertex
    xmax <- half_w * (1 - (ay - half_h / 2) / (half_h / 2))
  }
  c(-xmax, xmax)
}

# --- Generate subplot data: two series with dynamic correlation ---
set.seed(42)

n <- 120
tt <- seq(0, 4 * pi, length.out = n)

# Smoothly varying correlation
rho <- 0.3 * sin(tt * 0.8) + 0.2 * cos(tt * 0.4)

# Generate correlated series
z1 <- sin(tt * 1.2) + 0.3 * sin(tt * 2.5) + rnorm(n, 0, 0.15)
z2 <- rho * z1 + sqrt(1 - rho^2) * (sin(tt * 1.2 + 1.5) + 0.3 * cos(tt * 3) + rnorm(n, 0, 0.15))

# Smooth for visual appeal
z1s <- stats::filter(z1, rep(1 / 5, 5), sides = 2) |> as.numeric()
z2s <- stats::filter(z2, rep(1 / 5, 5), sides = 2) |> as.numeric()

# Rescale helper
rescale <- function(x, lo, hi) {
  r <- range(x, na.rm = TRUE)
  lo + (x - r[1]) / (r[2] - r[1]) * (hi - lo)
}

# Position the time series in the centre of the hex
x_vals  <- rescale(tt[!is.na(z1s)], -0.52, 0.52)
y1_vals <- rescale(z1s[!is.na(z1s)], -0.35, 0.35)
y2_vals <- rescale(z2s[!is.na(z2s)], -0.35, 0.35)

line_df <- data.frame(
  x = rep(x_vals, 2),
  y = c(y1_vals, y2_vals),
  series = rep(c("Y1", "Y2"), each = length(x_vals))
)

# --- Marginal density curves (rotated 90deg, on left and right edges) ---
d1 <- density(z1s[!is.na(z1s)], n = 256, adjust = 1.2)
d2 <- density(z2s[!is.na(z2s)], n = 256, adjust = 1.2)

# Map density x-axis (data values) into the same vertical range as the time series
d1_y <- rescale(d1$x, -0.35, 0.35)
d2_y <- rescale(d2$x, -0.35, 0.35)

# Map density heights to horizontal width; anchor at the series edges
margin_width <- 0.25
d1_x_base <- -0.56
d2_x_base <-  0.56

d1_x_raw <- d1_x_base - rescale(d1$y, 0, margin_width)
d2_x_raw <- d2_x_base + rescale(d2$y, 0, margin_width)

# Clip density polygons to the hex boundary
clip_to_hex <- function(x_density, y_density, x_base_val, side = "left",
                        radius = 0.88) {
  x_out <- x_density
  x_base_out <- rep(x_base_val, length(y_density))
  for (i in seq_along(y_density)) {
    bounds <- hex_x_at_y(y_density[i], radius = radius)
    if (side == "left") {
      x_out[i] <- max(x_out[i], bounds[1])
      x_base_out[i] <- max(x_base_out[i], bounds[1])
    } else {
      x_out[i] <- min(x_out[i], bounds[2])
      x_base_out[i] <- min(x_base_out[i], bounds[2])
    }
  }
  # Filter out points where density is fully clipped
  keep <- if (side == "left") x_out < x_base_out else x_out > x_base_out
  list(x_density = x_out[keep], x_base = x_base_out[keep], y = y_density[keep])
}

cl1 <- clip_to_hex(d1_x_raw, d1_y, d1_x_base, "left")
cl2 <- clip_to_hex(d2_x_raw, d2_y, d2_x_base, "right")

margin1_df <- data.frame(
  x = c(cl1$x_density, rev(cl1$x_base)),
  y = c(cl1$y, rev(cl1$y))
)
margin2_df <- data.frame(
  x = c(cl2$x_density, rev(cl2$x_base)),
  y = c(cl2$y, rev(cl2$y))
)

# --- Colours ---
col_fill   <- "#023047"
col_border <- "#219EBC"
col_line1  <- "#8ECAE6"
col_line2  <- "#FFB703"

# --- Build the sticker ---
hex <- hex_polygon(0, 0, radius = 1)
hex_inner <- hex_polygon(0, 0, radius = 0.97)

p <- ggplot() +
  # Hex fill
  geom_polygon(data = hex, aes(x, y), fill = col_fill, colour = NA) +
  # Marginal density — left (Y1)
  geom_polygon(data = margin1_df, aes(x, y),
               fill = NA, colour = col_line1, linewidth = 0.6) +
  # Marginal density — right (Y2)
  geom_polygon(data = margin2_df, aes(x, y),
               fill = NA, colour = col_line2, linewidth = 0.6) +
  # Time series lines
  geom_line(data = line_df, aes(x = x, y = y, colour = series),
            linewidth = 0.65, alpha = 0.9) +
  scale_colour_manual(values = c("Y1" = col_line1, "Y2" = col_line2)) +
  # Hex border (on top)
  geom_polygon(data = hex_inner, aes(x, y),
               fill = NA, colour = col_border, linewidth = 2) +
  # Clean up
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# --- Save ---
outfile <- file.path(here::here(), "man", "figures", "logo.png")
ggsave(outfile, p, width = 2, height = 2 * 2 / sqrt(3),
       dpi = 300, bg = "transparent")

message("Hex sticker saved to ", outfile)
