

#this loop does not work, some ID's have few observations. Need to remove them.
for (i in 1/40) {
  test <- subset(from_psych_weather, ID==i)
  
  testframe <- data.frame(time = test$measurenr,
                          painintensity=test$painintensity,
                          painunpleasantness=test$painunpleasantness,
                          pom=test$barometric_pres,
                          rh=test$humidity,
                          tam=test$Temperature)
  acf2(testframe$painintensity, max.lag = 1)
}



test <- subset(from_psych_weather, ID=="11")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="29")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="6")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)



acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="15")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="7")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="8")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)
acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)

test <- subset(from_psych_weather, ID=="10")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)

acf2(testframe$painintensity)
acf2(testframe$painunpleasantness)
test <- subset(from_psych_weather, ID=="7")

testframe <- data.frame(time = test$measurenr,
                        painintensity=test$painintensity,
                        painunpleasantness=test$painunpleasantness,
                        pom=test$barometric_pres,
                        rh=test$humidity,
                        tam=test$Temperature)

acf(testframe$painintensity)
acf2(testframe$painunpleasantness)