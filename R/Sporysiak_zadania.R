library(agricolae)
library(MASS)
library(ks)
library(ggplot2)
library(class)
library(caret)
library(cluster)
library(mclust)
library(factoextra)
library(mixtools)


data(yacon)
df <- yacon[, c("locality", "stalks", "wff")]

table(df$locality)

#wykres
plot(df$stalks, df$wff, col = as.factor(df$locality), pch = 19, 
     xlab = "Stalks", ylab = "WFF", main = "Wykres punktowy")
#nie pomijałem tego jednego podejrzanego punktu

kde_list <- list()
levels_loc <- unique(df$locality)

for (cls in levels_loc) {
  subset_data <- df[df$locality == cls, c("stalks", "wff")]
  kde_list[[cls]] <- kde2d(x = subset_data$stalks, y = subset_data$wff)
}
#obszary klasyfikacji
par(mfrow = c(1, length(levels_loc)))
for (cls in levels_loc) {
  contour(kde_list[[cls]], main = paste("Estymator jądrowy dla klasy", cls), xlab = "Stalks", ylab = "WFF")
}

kda_model <- kda(x = df[, c("stalks", "wff")], x.group = df$locality)

table(kda_model$x.group, kda_model$x.group.estimate)

correct_reclassification <- sum(diag(table(kda_model$x.group, kda_model$x.group.estimate))) / nrow(df)
print(correct_reclassification)
# 0.9699074 poprawnej klasyfikacji
#b)

#dopasowanie rozkładu normalnego dla każdej klasy
norm_params <- list()
for (cls in levels_loc) {
  subset_data <- df[df$locality == cls, c("stalks", "wff")]
  norm_params[[cls]] <- list(mean = colMeans(subset_data), cov = cov(subset_data))
}

#funkcja klasyfikacji na podstawie największego prawdopodobieństwa
dmvnorm_wrapper <- function(x, mean, cov) {
  dmvnorm(x, mean, cov)
}

predict_class <- function(x, y) {
  probs <- sapply(levels_loc, function(cls) {
    dmvnorm_wrapper(c(x, y), norm_params[[cls]]$mean, norm_params[[cls]]$cov)
  })
  return(levels_loc[which.max(probs)])
}

df$predicted <- mapply(predict_class, df$stalks, df$wff)

conf_matrix <- table(df$locality, df$predicted)
print(conf_matrix)

correct_reclassification_norm <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(correct_reclassification_norm)
# 0.9259259 poprawnej klasyfikacji

#wykres z obszarami klasyfikacji
ggplot(df, aes(x = stalks, y = wff, color = as.factor(predicted))) +
  geom_point() +
  labs(title = "Obszary klasyfikacji na podstawie rozkładu normalnego") +
  theme_minimal()

#podpunkt c

# 0.9699074 poprawnej klasyfikacji dla estymatora jądroej gęstości
# 0.9259259 poprawnej klasyfikacji dla gęstości rozkładu normalnego
#lepszy wynik uzyskaliśmuy za pomocą estymatora jądrowwej gęstości od gęstości rozkładu normalnego


# Podpunkt d

#siatka punktów
x_range <- seq(min(df$stalks), max(df$stalks), length.out = 100)
y_range <- seq(min(df$wff), max(df$wff), length.out = 100)
grid <- expand.grid(stalks = x_range, wff = y_range)

#klasyfikacja knn dla k=3
knn_pred <- knn(train = df[, c("stalks", "wff")], test = grid, cl = df$locality, k = 3)
grid$predicted <- knn_pred

#obszary klasyfikacji
ggplot() +
  geom_tile(data = grid, aes(x = stalks, y = wff, fill = as.factor(predicted)), alpha = 0.3) +
  geom_point(data = df, aes(x = stalks, y = wff, color = as.factor(locality))) +
  labs(title = "Obszary klasyfikacji dla k-NN (k=3)") +
  theme_minimal()



#knn + loocv oraz 10-krotna krosswalidacja
set.seed(555)
k_values <- 1:10
k_results_loocv <- list()
k_results_cv10 <- list()

for (k in k_values) {
#loocv
  train_control_loocv <- trainControl(method = "LOOCV")
  knn_model_loocv <- train(locality ~ stalks + wff, data = df, method = "knn", 
                           trControl = train_control_loocv, tuneGrid = data.frame(k = k))
  k_results_loocv[[as.character(k)]] <- knn_model_loocv$results$Accuracy
  
#10-krotna kroswalidacja
  train_control_cv10 <- trainControl(method = "cv", number = 10)
  knn_model_cv10 <- train(locality ~ stalks + wff, data = df, method = "knn", 
                          trControl = train_control_cv10, tuneGrid = data.frame(k = k))
  k_results_cv10[[as.character(k)]] <- knn_model_cv10$results$Accuracy
}

print("Wyniki LOOCV:")
print(k_results_loocv)
print("Wyniki 10-krotnej kroswalidacji:")
print(k_results_cv10)
#najlepszy wynik uzyskaliśmy dla k=1 a następnie dla k=5

#e)
#dane
df <- read.csv("C:/Users/polik/Desktop/dane22.csv")
df_selected <- df[, c("A", "B")]

#wykresy
plot(df_selected, main = "Wykres punktowy danych", pch = 19)

cor(df_selected$A, df_selected$B)

#wybór liczby klastrów
klasy <- numeric(7)
for (i in 2:7) {
  klasy[i] <- kmeans(df_selected, centers = i, nstart = 25)$betweenss
}

#wykres betweenss
plot(1:7, klasy, type = "b", pch = 19, main = "Optymalna liczba klastrów (betweenss)", 
     xlab = "Liczba klastrów", ylab = "Betweenss")

#k-means dla k=2 i 3
km2 <- kmeans(df_selected, centers = 2, nstart = 25)
km3 <- kmeans(df_selected, centers = 3, nstart = 25)

#wykresy
par(mfrow = c(1,2))
plot(df_selected, col = km2$cluster, pch = 19, main = "K-means: 2 klastry")
points(km2$centers, col = "red", pch = 8, cex = 2)
plot(df_selected, col = km3$cluster, pch = 19, main = "K-means: 3 klastry")
points(km3$centers, col = "red", pch = 8, cex = 2)
par(mfrow = c(1,1))  # Reset układu wykresów

#EM
set.seed(555)
df_matrix <- as.matrix(df_selected) + jitter(as.matrix(df_selected), amount = 1e-4)

EM2 <- mvnormalmixEM(df_matrix, k = 2)
EM3 <- mvnormalmixEM(df_matrix, k = 3)

#przypisanie klastrów na podstawie największego prawdopodobieństwa
df_selected$cluster_EM2 <- apply(EM2$posterior, 1, which.max)
df_selected$cluster_EM3 <- apply(EM3$posterior, 1, which.max)

#wykresy klastrów dla EM
par(mfrow = c(1,2))
plot(df_selected$A, df_selected$B, col = df_selected$cluster_EM2, pch = 19, main = "EM: 2 klastry")
plot(df_selected$A, df_selected$B, col = df_selected$cluster_EM3, pch = 19, main = "EM: 3 klastry")
par(mfrow = c(1,1))

# Obliczanie wskaźnika sylwetki dla K-Means (2 i 3 klastry)
sil_kmeans_2 <- silhouette(km2$cluster, dist(df_selected))
sil_kmeans_3 <- silhouette(km3$cluster, dist(df_selected))

mean_sil_kmeans_2 <- mean(sil_kmeans_2[, 3])
mean_sil_kmeans_3 <- mean(sil_kmeans_3[, 3])

# Obliczanie wskaźnika sylwetki dla EM (2 i 3 klastry)
sil_em_2 <- silhouette(df_selected$cluster_EM2, dist(df_selected))
sil_em_3 <- silhouette(df_selected$cluster_EM3, dist(df_selected))

mean_sil_em_2 <- mean(sil_em_2[, 3])
mean_sil_em_3 <- mean(sil_em_3[, 3])

# Tworzenie tabeli porównawczej
comparison <- data.frame(
  Metoda = c("K-Means (2)", "K-Means (3)", "EM (2)", "EM (3)"),
  Wskaźnik_Sylwetki = c(mean_sil_kmeans_2, mean_sil_kmeans_3, mean_sil_em_2, mean_sil_em_3)
)

print(comparison)

#dla k=2 lepiej wypada EM, natomiast lepsze dopasowanie jest dla 3 klastrów w obu przypadkach i jest ono identyczne, dane są dobrze oddzielone i wyniki mogą takie być


