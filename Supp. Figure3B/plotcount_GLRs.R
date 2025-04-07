
#Use plotCounts function to plot GLR2.9a, GLR2.9b and GLR3.1 expression level in different samples
data = plotCounts(dds, gene = "Niben101Scf01212g02011", intgroup = c("Short_name"), returnData = T)
#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))
ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR2.9a")

data = plotCounts(dds, gene = "Niben101Scf08670g00020", intgroup = c("Short_name"), returnData = T)
#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))
ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR2.9b")

data = plotCounts(dds, gene = "Niben101Scf00090g03004", intgroup = c("Short_name"), returnData = T)
#Turn your 'Short_name' column into a character vector
data$Short_name <- as.character(data$Short_name)
#Then turn it back into a factor with the levels in the correct order
data$Short_name <- factor(data$Short_name, levels=unique(data$Short_name))
ggplot(data, aes(x=interaction(Short_name), y=count, color=Short_name)) + 
  geom_jitter(size=2) + 
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("counts (normalized by size factor)")+
  ggtitle("NbGLR3.1")