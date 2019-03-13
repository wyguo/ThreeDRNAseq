#' Generate HCL colours, which are the same to \code{ggplot} default colour palette
#' @param n number of colours to generate.
#' @return a vector of \code{n} colours.
#' @export
gg.color.hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate distinct colours
#' @param n number of colours to generate.
#' @return a vector of \code{n} colours.
#' @export
distinct.color <- function(n){
  col.lib <- c("#F8766D","#00BFC4","#1874CD","red","yellow","green","darkgreen",
               gg.color.hue(10),
               "#C0D57E","#BF6CF4","#D435BB","#C2BBCB","#CAE678","#BEC3AD","#EC75F6","#F8DAD6","#D9FBD7","#46C3CF","#5FECDE","#AEE1F5","#AA67CC",
               "#81C1F1","#D6A381","#81DADB","#99AD6D","#B8ECF3","#828D6F","#F832A9","#8DD092","#7FC0D0","#6EDDD0","#D4BE95","#B8F7F2","#9752C4",
               "#AF9D7C","#468CE2","#9AA843","#4F9F7A","#8C7DB6","#CA44D8","#93D0B9","#4A718E","#F259C8","#E1C358","#BB924E","#5682F1","#B3E18E",
               "#59C1A6","#3BA7EB","#4ED06A","#A28FF0","#F4F2E7","#F7C2DA","#6921ED","#E8AAE6","#C2F2DF","#A3F7A1","#F0A85C","#BFF399","#E93B87",
               "#D7C7F0","#BCF7B6","#39B778","#DD714B","#B7D9A6","#E19892","#BF70C4","#B882F4","#74F38B","#D6E3C2","#9E23B2","#65F55F","#97F3EB",
               "#6BCD37","#8A7B8E","#DEEDA7","#F03E4F","#8FE77B","#ACF165","#8FCAC4","#CA1FF5","#E854EE","#43E11C","#6EEDA7","#B584B3","#BAAAF2",
               "#DE88C9","#415B96","#4EF4F3","#81F8CA","#C1B4AC","#BBA4CF","#AB55E5","#F789EC","#CEDFF8","#A9BEF0","#658DC6","#94ACEA","#91CD6F",
               "#BAA431","#8FB2C5","#F8F1A3","#D1CCCE","#F533D8","#844C56","#F3DA7D","#BDE5C0","#F1F69C","#F17C85","#C5E743","#B8C8ED","#3ED3EC",
               "#E6C5D9","#AD7967","#EAF4C3","#3F4DC1","#F4DDEC","#F1B136","#C247F3","#6FD3A1","#DCF7A3","#B8D1C2","#545CED","#2C7A29","#E4B1BB",
               "#785095","#8373F1","#D0F2C5","#B53757","#49F79F","#F0BE8B","#9630DD","#D6D89F","#79EED3","#EEEFD7","#7FA41F","#D1BF71","#6AEFF6",
               "#D661D8","#CFEEF8","#9ECB59","#67D382","#A8945E","#BB9B9E","#F3E4B9","#9964F3","#835BD4","#A5E2BF","#9CA6C8","#CD9BE9","#8FF038",
               "#B6C08A","#D19BB5","#8893EB","#51746E","#5920AC","#F3A1BB","#D7E1DA","#D668C3","#EC9EF2","#C2738A","#EDE979","#615AB4","#E9EF68",
               "#F3DECA","#90A6A0","#F2B8E7","#36DEA3","#EED89A","#5AADA5","#92EFAD","#A3288A","#ECC9F4","#91C9E7","#ECC2B8","#EB9BCC","#2735BB",
               "#5E7A2D","#C3BDF7","#DCF8F4","#52E0BC","#C2F8D6","#B3CFDA","#6E93A2","#9AECCA","#E3D8F0","#D97EE7","#9CBE9A","#60A858","#55F9D7",
               "#F47AA2","#DCF5E3","#B55398","#AC81D3","#644AF2","#546DEA","#E7EB37","#EA7E2A","#A1517E","#E131E4","#DE92F4","#EC66B1","#82DDF3",
               "#E99475","#E4EAF3","#E8D144","#40C0EE","#559ABA")
  col.lib <- unique(col.lib)
  if(n>length(col.lib))
    stop(paste0('n should be smaller than ',length(col.lib)))
  col.lib[1:n]
}

# gg_legend <- function(g){ 
#   tmp <- ggplot_gtable(ggplot_build(g)) 
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#   legend <- tmp$grobs[[leg]] 
#   return(legend)
#   } 
