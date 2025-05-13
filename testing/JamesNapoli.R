#2024-09-05 James Napoli 

devtools::load_all()

kuru <- list.files(DOWNLOADS, "Kuru*", full.names = TRUE)
chars <- ReadTNTCharacters(kuru)

dataset <- chars
indices <- 823:824

dc <- Decompose(chars, 823:824)
attributes(dc)$or
dataset <- dc
MaximizeParsimony(dc)
