library(drc)
data(methionine)
#-----------------------------------------------------------------------------------------------
# Zadanie 1+2)
head(methionine)
plot(methionine$dose, methionine$gain, pch = 19, col = "darkgray",main = "Wykres danych",xlab = "Dawka - dose", ylab = "Przyrost - gain")
legend("bottomright", legend = c("Model A - Asymptotyczny", "Model B - Logistyczny"), col = c("blue", "red"), lwd = 2)
#Dane są nieliniowe, tzn ich układ przypomina logarytm albo pierwiastek, nie można zastosować przybliżania funkcją liniową, ponieważ MSE byłby zbyt duży
# Model A
#b1 oszacowane z wykresu (około 1450), b2 około 1700, b3 testowane wieeele razy
m_a = nls(gain ~ b1 + (b2 - b1) * (1 - exp(-b3 * dose)), data = methionine, start = list(b1 = 1450, b2 = 1700, b3 = 15))
curve(predict(m_a, data.frame(dose = x)), add = TRUE, col = "blue", lwd = 2)
summary(m_a)
#nls poprawił wartości b1 na 1455.472, b2 na 1716.149, i b3 na 15.022. RSE=18.18

# Model B
#b3 oszacowane z wykresu, b1 i b2 analogicznie jak b3 powyżej
m_b = nls(gain ~ b3 / (1 + exp(b1 + b2 * dose)), data = methionine, start = list(b1 = -2, b2 = -15, b3 = 1700))
curve(predict(m_b, data.frame(dose = x)), add = TRUE, col = "red", lwd = 2)
summary(m_b)
#nls poprawił wartości b1 na -1.73621, b2 na -16.89654, i b3 na 1713.03248. RSE=18.44

n = nrow(methionine)
#MSE dla A
rss_a = sum(resid(m_a)^2)
mse_a = rss_a / n
cat("MSE dla A:", mse_a) #MSE dla A: 220.3054

#MSE dla B
rss_b = sum(resid(m_b)^2)
mse_b = rss_b / n
cat("MSE dla B:", mse_b) #MSE dla B: 226.6285

#W obu przypadkach i z MSE i RSE mniejszy błąd jest dla modelu A, więc jest on lepiej dopasowany dla istniejących danych
#Graficznie różnica istnieje, ale na oko jest ona niezauważalna pod względem 'która jest lepsza'. są po prostu bardzo podobne
#----------------------------------------------------------------------------------------------------------
#Zadanie 3
mod_loess = loess(gain ~ dose, data = methionine)

#RSS (nieparametryczne, więc nie będzie tu żadnych parametrów podwwanych)
RSS_loess = sum(mod_loess$residuals^2)
cat("RSS dla LOESS:", RSS_loess) #RSS dla LOESS: 1886.5
#model dużo gorzej poradził sobie z przybliżeniem funkcji, błąd jest bardzo duży w tym przypadku

#Spliny
mod_spline = smooth.spline(methionine$dose, methionine$gain)
#RSS między faktycznym, a wyliczonym przez spliny
y_pred_spline = predict(mod_spline, methionine$dose)$y
RSS_spline = sum((methionine$gain - y_pred_spline)^2)
cat("RSS dla Splinów:", RSS_spline) #RSS dla Splinów: 2348.383
#Model poradził sobie gorzej niż LOESS

plot(methionine$dose, methionine$gain, pch = 19, col = "grey",
     main = "LOESS i Spline", xlab = "Dawka", ylab = "Przyrost")
lines(methionine$dose, predict(mod_loess), col = "purple", lwd = 2)
lines(mod_spline, col = "darkgreen", lwd = 2)
legend("bottomright", legend = c("LOESS", "Spliny"), col = c("purple", "darkgreen"), lwd = 2)
#Loess stworzył 2 linie na wykresie, RSE policzone jest dobrze, ale źle narysowane, przez kolejność danych (linia złączyła slę dla jednego punktu jako wejście i wyjście, ważna jest ta górna).
#--------------------------------------------------------------------------------------------------------------
#Zadanie 4
modelA_init = function(mCall, LHS, data) {xy = sortedXyData(mCall[["dose"]], LHS, data)
  b1 = min(xy$y)
  b2 = max(xy$y)
  b3 = 10
  value = c(b1, b2, b3)
  names(value) = c("b1", "b2", "b3")
  return(value)
}
SSmodelA = selfStart(~ b1 + (b2 - b1) * (1 - exp(-b3 * dose)), initial = modelA_init, parameters = c("b1", "b2", "b3"))
model_test = nls(gain ~ SSmodelA(dose, b1, b2, b3), data = methionine)
summary(model_test)
#B1: 1455.472, B2: 1716.149, B3: 15.022, RSE: 18.18 
#Wartości są identyczne jak w modelu A, więc obie funkcje znalażły to samo 'optymalne' rowziązanie
#-----------------------------------------------------------------------------------------------------
# Zadanie 5
bledy = c()

#pętla n razy dla krosswalidacji n-krotnej
for(i in 1:n) {
  #LOOCV
  m = nls(gain ~ SSmodelA(dose, b1, b2, b3), data = methionine[-i,])
  #i-ty błąd ((y) rzecziwisty - przewidziany)^2
  bledy[i] = (methionine$gain[i] - predict(m, methionine[i,]))^2
}

RSS_CV = sum(bledy)
cat("Wynik kroswalidacji:", RSS_CV) #Wynik Kroswalidacji: 4494.492
#Błąd jest dużo większy, ponieważ usunięcie jednego punktu ze zbioru danych zawierającego ich 9 ma bardzo duży wpływ na
#wartości błędu. Największe znaczenie miał tam 1 punkt, bo był najdalej oddalony od pozostałyhc, wpływającmocno na wartości
------------------------------------------------------------------------------------------------------
#ogólnie zrobione zadania 1-5, bez 6 i 7
