
# library(cowplot)
# i1 <- ggdraw() + draw_image("figures/p400_c.c_c.w.png")
# i2 <- ggdraw() + draw_image("figures/p800_c.c_c.w.png")
# i3 <- ggdraw() + draw_image("figures/p400_w.w_w.c.png")
# i4 <- ggdraw() + draw_image("figures/p800_w.w_w.c.png")
# plot_grid(i1,i2,i3,i4, labels = "AUTO")
# cowplot::plot_grid(p400_cc_cw,p800_cc_cw, 
#                    p400_w.w_w.c, p800_w.w_w.c)


library(magick)
i1 <- image_read("figures/p400_c.c_w.w_.png")
i2 <- image_read("figures/p800_c.c_w.w_.png")
i3 <- image_read("figures/p400_c.c_c.w.png")
i4 <- image_read("figures/p800_c.c_c.w.png")
i5 <- image_read("figures/p400_w.w_w.c.png")
i6 <- image_read("figures/p800_w.w_w.c.png")
i1 <- image_annotate(i1, "(a)", size = 100, color = "black",
               location = "+190+100")
i2 <- image_annotate(i2, "(b)", size = 100, color = "black",
                     location = "+125+100")
i3 <- image_annotate(i3, "(c)", size = 100, color = "black",
                     location = "+190+100")
i4 <- image_annotate(i4, "(d)", size = 100, color = "black",
                     location = "+125+100")
i5 <- image_annotate(i5, "(e)", size = 100, color = "black",
                     location = "+190+100")
i6 <- image_annotate(i6, "(f)", size = 100, color = "black",
                     location = "+125+100")

img_top <- c(i1,i2)
img_mid <- c(i3,i4)
img_bottom <- c(i5,i6)
ijoin_top <- image_append(img_top)
ijoin_mid <- image_append(img_mid)
ijoin_bottom <- image_append(img_bottom)
ijoin_P400P800 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(ijoin_P400P800, path = "figures/JoinedP400P800_v2.png", format = "png")


library(magick)
i1 <- image_read("figures/vcmax400_c.c_w.w.png")
i2 <- image_read("figures/jmax400_c.c_w.w.png")
i3 <- image_read("figures/vcmax400_c.c_c.w.png")
i4 <- image_read("figures/jmax400_c.c_c.w.png")
i5 <- image_read("figures/vcmax400_w.w_w.c.png")
i6 <- image_read("figures/jmax400_w.w_w.c.png")
i1 <- image_annotate(i1, "(a)", size = 100, color = "black",
                     location = "+190+40")
i2 <- image_annotate(i2, "(b)", size = 100, color = "black",
                     location = "+190+40")
i3 <- image_annotate(i3, "(c)", size = 100, color = "black",
                     location = "+190+40")
i4 <- image_annotate(i4, "(d)", size = 100, color = "black",
                     location = "+190+40")
i5 <- image_annotate(i5, "(e)", size = 100, color = "black",
                     location = "+190+30")
i6 <- image_annotate(i6, "(f)", size = 100, color = "black",
                     location = "+190+30")

img_top <- c(i1,i2)
img_mid <- c(i3,i4)
img_bottom <- c(i5,i6)
ijoin_top <- image_append(img_top)
ijoin_mid <- image_append(img_mid)
ijoin_bottom <- image_append(img_bottom)
ijoin_P400P800 <- image_append(c(ijoin_top,ijoin_mid, ijoin_bottom), stack=T)

image_write(ijoin_P400P800, path = "figures/Joined_Vcmax_Jmax.png", format = "png")
