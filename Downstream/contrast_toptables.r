library(dplyr)
df1 <- data.frame(gene=c("Clu", "Lum", "Prg4"), cluster=1:3, description=c("haha", "hehe", "hoho"))
df2 <- data.frame(gene=c("ENSMUSG0000001", "ENSMUSG0000002", "ENSMUSG0000003"), cluster=1:3, description=c("hihi", "hoho", "hehe"))

merge(df1, df2, by="description", all.x=TRUE)

df1 %>% inner_join(df2, by="description")
