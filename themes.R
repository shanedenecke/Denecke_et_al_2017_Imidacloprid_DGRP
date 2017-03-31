
spread_theme <-theme(text=element_text(face="bold",family="serif"),axis.line.x =element_blank(), axis.ticks.x=element_blank(),panel.grid=element_blank(),
                    panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=20),strip.background=element_rect("grey95"),
                    axis.title=element_text(size=17),legend.position="none",axis.text.x=element_blank())

corr_theme <- theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),axis.ticks.x=element_line(),
                    panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=20),strip.background=element_rect("grey95"),
                    axis.title=element_text(size=17),legend.position="none")


cyp6g1.ind_theme <-  theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),axis.ticks.x=element_line(),
                           panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=20),strip.background=element_rect("grey95"),
                           axis.title=element_text(size=17))



dgrp.hplc_theme <- theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),axis.ticks.x=element_line(),
                        panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=20),strip.background=element_rect("grey95"),
                        axis.title=element_text(size=17),axis.text.x=element_text(angle=45,hjust=1,size=8))




g1KO.wiggle_theme <- theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),axis.ticks.x=element_line(),
                         panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=18),strip.background=element_rect("grey95"),
                         axis.title=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1),legend.position = "none")


rt_theme <- theme(text=element_text(face="bold",family="serif"),panel.grid=element_blank(),axis.ticks.x=element_line(),
                           panel.border=element_rect(colour="black",fill=NA),strip.text=element_text(size=18),strip.background=element_rect("grey95"),
                           axis.title=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1))
