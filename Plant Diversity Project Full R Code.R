
library("vegan")
library("dplyr")
library("ggplot2")
library("data.table")
library("ggmap")
library("maps")
library("stringr")
library("car")

## Functions and Metadata

setwd("C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data")

AstL <- read.csv("Asterales.csv")
EriL <- read.csv("Ericales.csv")
AstL1 <- subset(AstL, select = -c(occurrenceID, datasetKey, kingdom, phylum, class, countryCode, occurrenceStatus, publishingOrgKey, coordinatePrecision, elevation, elevationAccuracy, depth, depthAccuracy, eventDate, day, month, year, taxonKey, speciesKey, institutionCode, catalogNumber, identifiedBy, dateIdentified, license, rightsHolder, recordedBy, typeStatus, lastInterpreted, mediaType, issue, individualCount, basisOfRecord, collectionCode, recordNumber, establishmentMeans))
EriL1 <- subset(EriL, select = -c(occurrenceID, datasetKey, kingdom, phylum, class, countryCode, occurrenceStatus, publishingOrgKey, coordinatePrecision, elevation, elevationAccuracy, depth, depthAccuracy, eventDate, day, month, year, taxonKey, speciesKey, institutionCode, catalogNumber, identifiedBy, dateIdentified, license, rightsHolder, recordedBy, typeStatus, lastInterpreted, mediaType, issue, individualCount, basisOfRecord, collectionCode, recordNumber, establishmentMeans))

PF1 <- "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data\\GBIF1"
PF2 <- "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data\\GBIF2"
PF3 <- "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data\\GBIF3"
PF4 <- "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data\\GBIF4"
PF5 <- "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Raw_Data\\GBIF5"

combine_df <- function(x){
  csv_files <- list.files(x, pattern = "\\.csv$", full.names = TRUE)
  list_of_dataframes <- lapply(csv_files, read.csv)
  combined_data <- do.call(rbind, list_of_dataframes)
  combined_data_clean <- subset(combined_data, select = -c(occurrenceID, datasetKey, kingdom, phylum, class, countryCode, occurrenceStatus, publishingOrgKey, coordinatePrecision, elevation, elevationAccuracy, depth, depthAccuracy, eventDate, day, month, year, taxonKey, speciesKey, institutionCode, catalogNumber, identifiedBy, dateIdentified, license, rightsHolder, recordedBy, typeStatus, lastInterpreted, mediaType, issue, individualCount, basisOfRecord, collectionCode, recordNumber, establishmentMeans))
}

rdf1 <- combine_df(PF1)
rdf2 <- combine_df(PF2)
rdf3 <- combine_df(PF3)
rdf4 <- combine_df(PF4)
rdf5 <- combine_df(PF5)

create_df <- function(a, b){
  df <- rbind(AstL1[AstL1$decimalLatitude > a & AstL1$decimalLatitude < b, ], EriL1[EriL1$decimalLatitude > a & EriL1$decimalLatitude < b, ],  
              rdf1[rdf1$decimalLatitude > a & rdf1$decimalLatitude < b, ], 
              rdf2[rdf2$decimalLatitude > a & rdf2$decimalLatitude < b, ], 
              rdf3[rdf3$decimalLatitude > a & rdf3$decimalLatitude < b, ], 
              rdf4[rdf4$decimalLatitude > a & rdf4$decimalLatitude < b, ], 
              rdf5[rdf5$decimalLatitude > a & rdf5$decimalLatitude < b, ])
}
df32 <- create_df(29.999999, 32.000001)
df34 <- create_df(32.000000, 34.000001)
df36 <- create_df(34.000000, 36.000001)
df38 <- create_df(36.000000, 38.000001)
df40 <- create_df(38.000000, 40.000001)
df42 <- create_df(40.000000, 42.000001)
df44 <- create_df(42.000000, 44.000001)
df46 <- create_df(44.000000, 46.000001)
df48 <- create_df(46.000000, 48.000001)
df50 <- create_df(48.000000, 50.000001)

jfam <- function(x, y) { 
  return (n_distinct(intersect(x$family, y$family))/ 
            n_distinct(union(x$family, y$family))) 
}

jgen <- function(x, y) { 
  return (n_distinct(intersect(x$genus, y$genus))/ 
            n_distinct(union(x$genus, y$genus))) 
}

jspec <- function(x, y) { 
  return (n_distinct(intersect(x$species, y$species))/ 
            n_distinct(union(x$species, y$species))) 
}

create_dflat2 <- function(a, b, c, d, e, f, g, h, i) {
  return(dflat2 <- data.frame(w = c(log10(jfam(a, b)), log10(jfam(b, c)), log10(jfam(c, d)), log10(jfam(d, e)), log10(jfam(e, f)), log10(jfam(f, g)), log10(jfam(g, h)), log10(jfam(h, i))),
                              x = c(log10(jgen(a, b)), log10(jgen(b, c)), log10(jgen(c, d)), log10(jgen(d, e)), log10(jgen(e, f)), log10(jgen(f, g)), log10(jgen(g, h)), log10(jgen(h, i))), 
                              y = c(log10(jspec(a, b)), log10(jspec(b, c)), log10(jspec(c, d)), log10(jspec(d, e)), log10(jspec(e, f)), log10(jspec(f, g)), log10(jspec(g, h)), log10(jspec(h, i))),
                              z = c(34, 36, 38, 40, 42, 44, 46, 48)))
}
dflat2 <- create_dflat2(df32, df34, df36, df38, df40, df42, df44, df46, df48)
dflat2
LatLongdf <- function (b, c, d, e, f){
  return (data.frame(b[b$decimalLatitude >= e & b$decimalLatitude <= c & b$decimalLongitude <= f & b$decimalLongitude >= d, ]))
}  


graph_gen <- ggplot(dflat2, aes(x=z, y=x)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  annotate("text",x=36.5,y=-0.05,label=(paste0("y=0.002529636x-0.234527874"))) +
  xlab("Comparison of diversity across n degrees latitude") + 
  ylab("ln of variation of genera between localities") +
  ylim(-0.3, -0.05) +
  ggtitle("Variation in Genera in North America Accross on a Latitudinal Gradient")
graph_gen

graph_spec <- ggplot(dflat2, aes(x=z, y=y)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  annotate("text",x=36.5,y=-0.05,label=(paste0("y=0.004624983x-0.418402136"))) +
  xlab("Comparison of diversity across n degrees latitude") + 
  ylab("ln of variation of species between localities") +
  ylim(-0.3, -0.05) +
  ggtitle("Variation in Species in North America Accross on a Latitudinal Gradient")
graph_spec

graph_fam <- ggplot(dflat2, aes(x=z, y=w)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  annotate("text",x=36.5,y=0,label=(paste0("y=-0.000649038x-0.037930638"))) +
  xlab("Comparison of diversity across n degrees latitude") + 
  ylab("ln of variation of families between localities") +
  ylim(-0.2, -0) +
  ggtitle("Variation in Families in North America Accross on a Latitudinal Gradient")
graph_fam

lat_greater <- function(x, y, z) {
  return(x[x$decimalLatitude > y & x$stateProvince == z, ])
}

lat_less <- function(x, y, z) {
  return(x[x$decimalLatitude < y & x$stateProvince == z, ])
}

long_less <- function(x, y, z) {
  return(x[x$decimalLongitude < y & x$stateProvince == z, ])
}

long_greater <- function(x, y, z) {
  return(x[x$decimalLongitude > y & x$stateProvince == z, ])
}
undup <- function(a) {
  a[!duplicated(a$gbifID), ]
}

plot_data <- function(a){
  us_map <- map_data("state")
  p <- ggplot() +
    geom_polygon(data = us_map, aes(x = long, y = lat, group = group), 
                 fill = "white", color = "black", size = 0.5) +
    coord_fixed(ratio = 1.2, xlim = c(-125, -66), ylim = c(24, 49)) +
    theme_void()
  return(p + geom_point(data = a, aes(x = a$decimalLongitude, y = a$decimalLatitude), color = "red", size = 0.5))
}

plot_data2 <- function(a){
  us_map <- map_data("state")
  p <- ggplot() +
    geom_polygon(data = us_map, aes(x = long, y = lat, color = biome), 
                 fill = "white", color = "black", size = 0.5) +
    coord_fixed(ratio = 1.2, xlim = c(-125, -66), ylim = c(24, 49)) +
    theme_void()
  return(p + geom_point(data = a, aes(x = a$decimalLongitude, y = a$decimalLatitude), size = 0.25))
}

less_mlatlongdf <- function(a, b, c, d, e, f){
  x1 <- c
  y1 <- b
  x2 <- e
  y2 <- d
  m <- ((y2 - y1) / (x2 - x1))
  int <- y1 - (m * x1)
  subset <- (data.frame(a[a$decimalLatitude > b & a$decimalLatitude < d & a$decimalLatitude < ((m * a$decimalLongitude) + int) & a$decimalLongitude > ((a$decimalLatitude - int)/m), ]))
  return(subset[subset$stateProvince == f, ])
}
greater_mlatlongdf <- function(a, b, c, d, e, f){
  x1 <- c
  y1 <- b
  x2 <- e
  y2 <- d
  m <- ((y2 - y1) / (x2 - x1))
  int <- y1 - (m * x1)
  subset <- (data.frame(a[a$decimalLatitude > b &a$decimalLatitude <d & a$decimalLatitude > ((m * a$decimalLongitude) + int) & a$decimalLongitude < ((a$decimalLatitude - int)/m), ]))
  return(subset[subset$stateProvince == f, ])
}

ngreater_mlatlongdf <- function(a, b, c, d, e, f){
  x1 <- c
  y1 <- b
  x2 <- e
  y2 <- d
  m <- ((y2 - y1) / (x2 - x1))
  int <- y1 - (m * x1)
  subset <- (data.frame(a[a$decimalLatitude > d & a$decimalLatitude < b & a$decimalLatitude > ((m * a$decimalLongitude) + int) & a$decimalLongitude > ((a$decimalLatitude - int)/m), ]))
  return(subset[subset$stateProvince == f, ])
}
nless_mlatlongdf <- function(a, b, c, d, e, f){
  x1 <- c
  y1 <- b
  x2 <- e
  y2 <- d
  m <- ((y2 - y1) / (x2 - x1))
  int <- y1 - (m * x1)
  subset <- (data.frame(a[a$decimalLatitude > d & a$decimalLatitude < b & a$decimalLatitude < ((m * a$decimalLongitude) + int) & a$decimalLongitude < ((a$decimalLatitude - int)/m), ]))
  return(subset[subset$stateProvince == f, ])
}
draw_complex <- function(a, b, c, d, x, e, f, g, h, y, i, j, k, l, z, s){
  if(a-c > 0){
    if(a > 48 & c > 48) {
      r <- df50
    }else{ 
      if(a > 48 & c > 46 & c < 48){
        r <- df48_50
      }else{
        if(a < 48 & a > 46 & c > 46 & c < 48){
          r <- df48
        }else{
          if(a > 46 & a <48 & c < 46 & c > 44){
            r <- df46_48
          }else{
            if(a < 46 & a > 44 & c < 46 & c > 44){
              r<-df46
            }else{
              if(a < 46 & a > 44 & c < 44 & c > 42){
                r<-df44_46
              }else {
                if(a < 44 & a > 42 & c < 44 & c > 42){
                  r<-df44
                }else{
                  if(a < 44 & a > 42 & c < 42 & c > 40){
                    r <- df42_44
                  }else{
                    if(a < 42 & a > 40 & c < 42 & c > 40){
                      r<- df42
                    }else{
                      if(a < 42 & a > 40 & c < 40 & c > 38){
                        r<-df40_42
                      }else{
                        if(a < 40 & a > 38 & c < 40 & c > 38){
                          r<-df40
                        }else{
                          if(a < 40 & a > 38 & c < 38 & c > 36){
                            r <- df38_40
                          }else{
                            if(a < 38 & a > 36 & c < 38 & c > 36){
                              r <- df38
                            }else{
                              if(a < 38 & a > 36 & c < 36 & c > 34){
                                r <- df36_38
                              }else{ 
                                if(a < 36 & a > 34 & c < 36 & c > 34){
                                  r <- df36
                                }else{
                                  if(a < 36 & a > 34 & c < 34 & c > 32){
                                    r <- df34_36
                                  }else{
                                    if(a < 34 & a > 32 & c < 34 & c > 32){
                                      r<-df34
                                    }else{
                                      if(a < 34 & a > 32 & c < 32){
                                        r<-df32_34
                                      }else{
                                        r<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }else{
    if(c > 48 & a > 48) {
      r <- df50
    }else{ 
      if(c > 48 & a > 46 & a < 48){
        r <- df48_50
      }else{
        if(c < 48 & c > 46 & a > 46 & a < 48){
          r <- df48
        }else{
          if(c > 46 & c <48 & a < 46 & a > 44){
            r <- df46_48
          }else{
            if(c < 46 & c > 44 & a < 46 & a > 44){
              r<-df46
            }else{
              if(c < 46 & c > 44 & a < 44 & a > 42){
                r<-df44_46
              }else {
                if(c < 44 & c > 42 & a < 44 & a > 42){
                  r<-df44
                }else{
                  if(c < 44 & c > 42 & a < 42 & a > 40){
                    r <- df42_44
                  }else{
                    if(c < 42 & c > 40 & a < 42 & a > 40){
                      r<- df42
                    }else{
                      if(c < 42 & c > 40 & a < 40 & a > 38){
                        r<-df40_42
                      }else{
                        if(c < 40 & c > 38 & a < 40 & a > 38){
                          r<-df40
                        }else{
                          if(c < 40 & c > 38 & a < 38 & a > 36){
                            r <- df38_40
                          }else{
                            if(c < 38 & c > 36 & a < 38 & a > 36){
                              r <- df38
                            }else{
                              if(c < 38 & c > 36 & a < 36 & a > 34){
                                r <- df36_38
                              }else{ 
                                if(c < 36 & c > 34 & a < 36 & a > 34){
                                  r <- df36
                                }else{
                                  if(c < 36 & c > 34 & a < 34 & a > 32){
                                    r <- df34_36
                                  }else{
                                    if(c < 34 & c > 32 & a < 34 & a > 32){
                                      r<-df34
                                    }else{
                                      if(c < 34 & c > 32 & a < 32){
                                        r<-df32_34
                                      }else{
                                        r<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(e-g > 0){
    if(e > 48 & g > 48) {
      r2 <- df50
    }else{ 
      if(e > 48 & g > 46 & g < 48){
        r2 <- df48_50
      }else{
        if(e < 48 & e > 46 & g > 46 & g < 48){
          r2 <- df48
        }else{
          if(e > 46 & e <48 & g < 46 & g > 44){
            r2 <- df46_48
          }else{
            if(e < 46 & e > 44 & g < 46 & g > 44){
              r2<-df46
            }else{
              if(e < 46 & e > 44 & g < 44 & g > 42){
                r2<-df44_46
              }else {
                if(e < 44 & e > 42 & g < 44 & g > 42){
                  r2<-df44
                }else{
                  if(e < 44 & e > 42 & g < 42 & g > 40){
                    r2 <- df42_44
                  }else{
                    if(e < 42 & e > 40 & g < 42 & g > 40){
                      r2<- df42
                    }else{
                      if(e < 42 & e > 40 & g < 40 & g > 38){
                        r2<-df40_42
                      }else{
                        if(e < 40 & e > 38 & g < 40 & g > 38){
                          r2 <-df40
                        }else{
                          if(e < 40 & e > 38 & g < 38 & g > 36){
                            r2 <- df38_40
                          }else{
                            if(e < 38 & e > 36 & g < 38 & g > 36){
                              r2 <- df38
                            }else{
                              if(e < 38 & e > 36 & g < 36 & g > 34){
                                r2 <- df36_38
                              }else{ 
                                if(e < 36 & e > 34 & g < 36 & g > 34){
                                  r2 <- df36
                                }else{
                                  if(e < 36 & e > 34 & g < 34 & g > 32){
                                    r2 <- df34_36
                                  }else{
                                    if(e < 34 & e > 32 & g < 34 & g > 32){
                                      r2<-df34
                                    }else{
                                      if(e < 34 & e > 32 & g < 32){
                                        r2<-df32_34
                                      }else{
                                        r2<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }else{
    if(g > 48 & e > 48) {
      r2 <- df50
    }else{ 
      if(g > 48 & g > 46 & e < 48){
        r2 <- df48_50
      }else{
        if(g < 48 & g > 46 & e > 46 & e < 48){
          r2 <- df48
        }else{
          if(g > 46 & g <48 & e < 46 & e > 44){
            r2 <- df46_48
          }else{
            if(g < 46 & g > 44 & e < 46 & e > 44){
              r2<-df46
            }else{
              if(g < 46 & g > 44 & e < 44 & e > 42){
                r2<-df44_46
              }else{
                if(g < 44 & g > 42 & e < 44 & e > 42){
                  r2<-df44
                }else{
                  if(g < 44 & g > 42 & e < 42 & e > 40){
                    r2 <- df42_44
                  }else{
                    if(g < 42 & g > 40 & e < 42 & e > 40){
                      r2<- df42
                    }else{
                      if(g < 42 & g > 40 & e < 40 & e > 38){
                        r2<-df40_42
                      }else{
                        if(g < 40 & g > 38 & e < 40 & e > 38){
                          r2 <-df40
                        }else{
                          if(g < 40 & g > 38 & e < 38 & e > 36){
                            r2 <- df38_40
                          }else{
                            if(g < 38 & g > 36 & e < 38 & e > 36){
                              r2 <- df38
                            }else{
                              if(g < 38 & g > 36 & e < 36 & e > 34){
                                r2 <- df36_38
                              }else{ 
                                if(g < 36 & g > 34 & e < 36 & e > 34){
                                  r2 <- df36
                                }else{
                                  if(g < 36 & g > 34 & e < 34 & e > 32){
                                    r2 <- df34_36
                                  }else{
                                    if(g < 34 & g > 32 & e < 34 & e > 32){
                                      r2<-df34
                                    }else{
                                      if(g < 34 & g > 32 & e < 32){
                                        r2<-df32_34
                                      }else{
                                        r2<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(i-k > 0){
    if(i> 48 & k > 48){
      r3 <- df50
    }else{ 
      if(i > 48 & k > 46 & k < 48){
        r3 <- df48_50
      }else{
        if(i < 48 & i > 46 & k > 46 & k < 48){
          r3 <- df48
        }else{
          if(i > 46 & i <48 & k < 46 & k > 44){
            r3 <- df46_48
          }else{
            if(i < 46 & i > 44 & k < 46 & k > 44){
              r3<-df46
            }else{
              if(i < 46 & i > 44 & k < 44 & k > 42){
                r3<-df44_46
              }else {
                if(i < 44 & i > 42 & k < 44 & k > 42){
                  r3<-df44
                }else{
                  if(i < 44 & i > 42 & k < 42 & k > 40){
                    r3 <- df42_44
                  }else{
                    if(i < 42 & i > 40 & k < 42 & k > 40){
                      r3<- df42
                    }else{
                      if(i < 42 & i > 40 & k < 40 & k > 38){
                        r3<-df40_42
                      }else{
                        if(i < 40 & i > 38 & k < 40 & k > 38){
                          r3<-df40
                        }else{
                          if(i < 40 & i > 38 & k < 38 & k > 36){
                            r3 <- df38_40
                          }else{
                            if(i < 38 & i > 36 & k < 38 & k > 36){
                              r3 <- df38
                            }else{
                              if(i < 38 & i > 36 & k < 36 & k > 34){
                                r3 <- df36_38
                              }else{ 
                                if(i < 36 & i > 34 & k < 36 & k > 34){
                                  r3 <- df36
                                }else{
                                  if(i < 36 & i > 34 & k < 34 & k > 32){
                                    r3 <- df34_36
                                  }else{
                                    if(i < 34 & i > 32 & k < 34 & k > 32){
                                      r3<-df34
                                    }else{
                                      if(i < 34 & i > 32 & k < 32){
                                        r3<-df32_34
                                      }else{
                                        r3<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }else{
    if(k > 48 & i > 48) {
      r3 <- df50
    }else{ 
      if(k > 48 & i > 46 & i < 48){
        r3 <- df48_50
      }else{
        if(k < 48 & k > 46 & i > 46 & i < 48){
          r3 <- df48
        }else{
          if(k > 46 & k <48 & i < 46 & i > 44){
            r3 <- df46_48
          }else{
            if(k < 46 & k > 44 & i < 46 & i > 44){
              r3<-df46
            }else{
              if(k < 46 & k > 44 & i < 44 & i > 42){
                r3<-df44_46
              }else {
                if(k < 44 & k > 42 & i < 44 & i > 42){
                  r3<-df44
                }else{
                  if(k < 44 & k > 42 & i < 42 & i > 40){
                    r3 <- df42_44
                  }else{
                    if(k < 42 & k > 40 & i < 42 & i > 40){
                      r3<- df42
                    }else{
                      if(k < 42 & k > 40 & i < 40 & i > 38){
                        r3<-df40_42
                      }else{
                        if(k < 40 & k > 38 & i  < 40 & i > 38){
                          r3<-df40
                        }else{
                          if(k < 40 & k > 38 & i < 38 & i > 36){
                            r3 <- df38_40
                          }else{
                            if(k < 38 & k > 36 & i < 38 & i > 36){
                              r3 <- df38
                            }else{
                              if(k < 38 & k > 36 & i < 36 & i > 34){
                                r3 <- df36_38
                              }else{ 
                                if(k < 36 & k > 34 & i < 36 & i > 34){
                                  r3 <- df36
                                }else{
                                  if(k < 36 & k > 34 & i  < 34 & i > 32){
                                    r3 <- df34_36
                                  }else{
                                    if(k < 34 & k > 32 & i < 34 & i > 32){
                                      r3<-df34
                                    }else{
                                      if(k < 34 & k > 32 & i < 32){
                                        r3<-df32_34
                                      }else{
                                        r3<-df32
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(x == "greater"){
    r4 <- greater_mlatlongdf(r, c, d, a, b, s)
  }else{
    if(x == "less"){
      r4 <- less_mlatlongdf(r, c, d, a, b, s)
    }else{
      if(x=="ngreater"){
        r4<-ngreater_mlatlongdf(r, c, d, a, b, s)
      }else{
        if(x == "nless"){
          r4<- nless_mlatlongdf(r, c, d, a, b, s)
        }
      }
    }
  }
  if(y == "greater"){
    r5<-greater_mlatlongdf(r2, g, h, e, f, s)
  }else{
    if(y == "less"){
      r5<-less_mlatlongdf(r2, g, h, e, f, s)
    }else{
      if(y == "ngreater"){
        r5<-ngreater_mlatlongdf(r2, g, h, e, f, s)
      }else{
        r5<- nless_mlatlongdf(r2, g, h, e, f, s)
      }
    }
  }
  if(z == "greater"){
    r6<-greater_mlatlongdf(r3, k, l, i, j, s)
  }else{
    if(z == "less"){
      r6<-less_mlatlongdf(r3, k, l, i, j, s)
    }else{
      if(z=="ngreater"){
        r6<-ngreater_mlatlongdf(r3, k, l, i, j, s)
      }else{
        r6<- nless_mlatlongdf(r3, k, l, i, j, s)
      }
    }
  }
  return(rbind(r4, r5, r6))
}

count_species <- function(a){
  nodup <- undup(a)
  species <- table(nodup$species)
  species_counts_df <- as.data.frame(species)
  return(species_counts_df)
}
gen_shannon <- function(a){
  r30.5 <- a[a$decimalLatitude < 30.5, ]
  cdf30.5 <- count_species(r30.5)
  if(nrow(cdf30.5) > 3){
    shannon_index30.5 <- diversity(cdf30.5[,-1], index = "shannon")
  }else{
    shannon_index30.5 <- c("N/A")
  }
  r31 <- a[a$decimalLatitude >= 30.5 & a$decimalLatitude < 31, ]
  cdf31 <- count_species(r31)
  if(nrow(cdf31) > 3){
    shannon_index31<- diversity(cdf31[,-1], index = "shannon")
  }else{
    shannon_index31 <- c("N/A")
  }
  r31.5 <- a[a$decimalLatitude >= 31 & a$decimalLatitude < 31.5, ]
  cdf31.5 <- count_species(r31.5)
  if(nrow(cdf31.5) > 3){
    shannon_index31.5 <- diversity(cdf31.5[,-1], index = "shannon")
  }else{
    shannon_index31.5<- c("N/A")
  }
  r32 <- a[a$decimalLatitude >= 31.5 & a$decimalLatitude < 32, ]
  cdf32 <- count_species(r32)
  if(nrow(cdf32) > 3){
    shannon_index32 <- diversity(cdf32[,-1], index = "shannon")
  }else{
    shannon_index32<- c("N/A")
  }
  
  r32.5 <- a[a$decimalLatitude >= 32 & a$decimalLatitude < 32.5, ]
  cdf32.5 <- count_species(r32.5)
  if(nrow(cdf32.5) > 3){
    shannon_index32.5<- diversity(cdf32.5[,-1], index = "shannon")
  }else{
    shannon_index32.5 <- c("N/A")
  }
  r33 <- a[a$decimalLatitude >= 32.5 & a$decimalLatitude < 33, ]
  cdf33 <- count_species(r33)
  if(nrow(cdf33) > 3){
    shannon_index33<- diversity(cdf33[,-1], index = "shannon")
  }else{
    shannon_index33 <- c("N/A")
  }
  r33.5 <- a[a$decimalLatitude >= 33 & a$decimalLatitude < 33.5, ]
  cdf33.5 <- count_species(r33.5)
  if(nrow(cdf33.5) > 3){
    shannon_index33.5<- diversity(cdf33.5[,-1], index = "shannon")
  }else{
    shannon_index33.5 <- c("N/A")
  }
  r34 <- a[a$decimalLatitude >= 33.5 & a$decimalLatitude < 34, ]
  cdf34 <- count_species(r34)
  if(nrow(cdf34) > 3){
    shannon_index34<- diversity(cdf34[,-1], index = "shannon")
  }else{
    shannon_index34 <- c("N/A")
  }
  r34.5 <- a[a$decimalLatitude >= 34 & a$decimalLatitude < 34.5, ]
  cdf34.5 <- count_species(r34.5)
  if(nrow(cdf34.5) > 3){
    shannon_index34.5<- diversity(cdf34.5[,-1], index = "shannon")
  }else{
    shannon_index34.5 <- c("N/A")
  }
  r35 <- a[a$decimalLatitude >= 34.5 & a$decimalLatitude < 35, ]
  cdf35 <- count_species(r35)
  if(nrow(cdf35) > 3){
    shannon_index35<- diversity(cdf35[,-1], index = "shannon")
  }else{
    shannon_index35 <- c("N/A")
  }
  r35.5 <- a[a$decimalLatitude >= 35 & a$decimalLatitude < 35.5, ]
  cdf35.5 <- count_species(r35.5)
  if(nrow(cdf35.5) > 3){
    shannon_index35.5<- diversity(cdf35.5[,-1], index = "shannon")
  }else{
    shannon_index35.5 <- c("N/A")
  }
  r36 <- a[a$decimalLatitude >= 35.5 & a$decimalLatitude < 36, ]
  cdf36 <- count_species(r36)
  if(nrow(cdf36) > 3){
    shannon_index36<- diversity(cdf36[,-1], index = "shannon")
  }else{
    shannon_index36 <- c("N/A")
  }
  r36.5 <- a[a$decimalLatitude >= 36 & a$decimalLatitude < 36.5, ]
  cdf36.5 <- count_species(r36.5)
  if(nrow(cdf36.5) > 3){
    shannon_index36.5<- diversity(cdf36.5[,-1], index = "shannon")
  }else{
    shannon_index36.5 <- c("N/A")
  }
  r37 <- a[a$decimalLatitude >= 36.5 & a$decimalLatitude < 37, ]
  cdf37 <- count_species(r37)
  if(nrow(cdf37) > 3){
    shannon_index37<- diversity(cdf37[,-1], index = "shannon")
  }else{
    shannon_index37 <- c("N/A")
  }
  r37.5 <- a[a$decimalLatitude >= 37 & a$decimalLatitude < 37.5, ]
  cdf37.5 <- count_species(r37.5)
  if(nrow(cdf37.5) > 3){
    shannon_index37.5<- diversity(cdf37.5[,-1], index = "shannon")
  }else{
    shannon_index37.5 <- c("N/A")
  }
  r38 <- a[a$decimalLatitude >= 37.5 & a$decimalLatitude < 38, ]
  cdf38 <- count_species(r38)
  if(nrow(cdf38) > 3){
    shannon_index38<- diversity(cdf38[,-1], index = "shannon")
  }else{
    shannon_index38 <- c("N/A")
  }
  r38.5 <- a[a$decimalLatitude >= 38 & a$decimalLatitude < 38.5, ]
  cdf38.5 <- count_species(r38.5)
  if(nrow(cdf38.5) > 3){
    shannon_index38.5<- diversity(cdf38.5[,-1], index = "shannon")
  }else{
    shannon_index38.5 <- c("N/A")
  }
  r39 <- a[a$decimalLatitude >= 38.5 & a$decimalLatitude < 39, ]
  cdf39 <- count_species(r39)
  if(nrow(cdf39) > 3){
    shannon_index39<- diversity(cdf39[,-1], index = "shannon")
  }else{
    shannon_index39 <- c("N/A")
  }
  r39.5 <- a[a$decimalLatitude >= 39 & a$decimalLatitude < 39.5, ]
  cdf39.5 <- count_species(r39.5)
  if(nrow(cdf39.5) > 3){
    shannon_index39.5<- diversity(cdf39.5[,-1], index = "shannon")
  }else{
    shannon_index39.5 <- c("N/A")
  }
  r40 <- a[a$decimalLatitude >= 39.5 & a$decimalLatitude < 40, ]
  cdf40 <- count_species(r40)
  if(nrow(cdf40) > 3){
    shannon_index40<- diversity(cdf40[,-1], index = "shannon")
  }else{
    shannon_index40 <- c("N/A")
  }
  r40.5 <- a[a$decimalLatitude >= 40 & a$decimalLatitude < 40.5, ]
  cdf40.5 <- count_species(r40.5)
  if(nrow(cdf40.5) > 3){
    shannon_index40.5<- diversity(cdf40.5[,-1], index = "shannon")
  }else{
    shannon_index40.5 <- c("N/A")
  }
  r41 <- a[a$decimalLatitude >= 40.5 & a$decimalLatitude < 41, ]
  cdf41 <- count_species(r41)
  if(nrow(cdf41) > 3){
    shannon_index41<- diversity(cdf41[,-1], index = "shannon")
  }else{
    shannon_index41 <- c("N/A")
  }
  r41.5 <- a[a$decimalLatitude >= 41 & a$decimalLatitude < 41.5, ]
  cdf41.5 <- count_species(r41.5)
  if(nrow(cdf41.5) > 3){
    shannon_index41.5<- diversity(cdf41.5[,-1], index = "shannon")
  }else{
    shannon_index41.5 <- c("N/A")
  }
  r42 <- a[a$decimalLatitude >= 41.5 & a$decimalLatitude < 42, ]
  cdf42 <- count_species(r42)
  if(nrow(cdf42) > 3){
    shannon_index42<- diversity(cdf42[,-1], index = "shannon")
  }else{
    shannon_index42 <- c("N/A")
  }
  r42.5 <- a[a$decimalLatitude >= 42 & a$decimalLatitude < 42.5, ]
  cdf42.5 <- count_species(r42.5)
  if(nrow(cdf42.5) > 3){
    shannon_index42.5<- diversity(cdf42.5[,-1], index = "shannon")
  }else{
    shannon_index42.5 <- c("N/A")
  }
  r43 <- a[a$decimalLatitude >= 42.5 & a$decimalLatitude < 43, ]
  cdf43 <- count_species(r43)
  if(nrow(cdf43) > 3){
    shannon_index43<- diversity(cdf43[,-1], index = "shannon")
  }else{
    shannon_index43 <- c("N/A")
  }
  r43.5 <- a[a$decimalLatitude >= 43 & a$decimalLatitude < 43.5, ]
  cdf43.5 <- count_species(r43.5)
  if(nrow(cdf43.5) > 3){
    shannon_index43.5<- diversity(cdf43.5[,-1], index = "shannon")
  }else{
    shannon_index43.5 <- c("N/A")
  }
  r44 <- a[a$decimalLatitude >= 43.5 & a$decimalLatitude < 44, ]
  cdf44 <- count_species(r44)
  if(nrow(cdf44) > 3){
    shannon_index44<- diversity(cdf44[,-1], index = "shannon")
  }else{
    shannon_index44 <- c("N/A")
  }
  r44.5 <- a[a$decimalLatitude >= 44 & a$decimalLatitude < 44.5, ]
  cdf44.5 <- count_species(r44.5)
  if(nrow(cdf44.5) > 3){
    shannon_index44.5<- diversity(cdf44.5[,-1], index = "shannon")
  }else{
    shannon_index44.5 <- c("N/A")
  }
  r45 <- a[a$decimalLatitude >= 44.5 & a$decimalLatitude < 45, ]
  cdf45 <- count_species(r45)
  if(nrow(cdf45) > 3){
    shannon_index45<- diversity(cdf45[,-1], index = "shannon")
  }else{
    shannon_index45 <- c("N/A")
  }
  r45.5 <- a[a$decimalLatitude >= 45 & a$decimalLatitude < 45.5, ]
  cdf45.5 <- count_species(r45.5)
  if(nrow(cdf45.5) > 3){
    shannon_index45.5<- diversity(cdf45.5[,-1], index = "shannon")
  }else{
    shannon_index45.5 <- c("N/A")
  }
  r46 <- a[a$decimalLatitude >= 45.5 & a$decimalLatitude < 46, ]
  cdf46 <- count_species(r46)
  if(nrow(cdf46) > 3){
    shannon_index46<- diversity(cdf46[,-1], index = "shannon")
  }else{
    shannon_index46 <- c("N/A")
  }
  r46.5 <- a[a$decimalLatitude >= 46 & a$decimalLatitude < 46.5, ]
  cdf46.5 <- count_species(r46.5)
  if(nrow(cdf46.5) > 3){
    shannon_index46.5<- diversity(cdf46.5[,-1], index = "shannon")
  }else{
    shannon_index46.5 <- c("N/A")
  }
  r47 <- a[a$decimalLatitude >= 46.5 & a$decimalLatitude < 47, ]
  cdf47 <- count_species(r47)
  if(nrow(cdf47) > 3){
    shannon_index47<- diversity(cdf47[,-1], index = "shannon")
  }else{
    shannon_index47 <- c("N/A")
  }
  r47.5 <- a[a$decimalLatitude >= 47 & a$decimalLatitude < 47.5, ]
  cdf47.5 <- count_species(r47.5)
  if(nrow(cdf47.5) > 3){
    shannon_index47.5<- diversity(cdf47.5[,-1], index = "shannon")
  }else{
    shannon_index47.5 <- c("N/A")
  }
  r48 <- a[a$decimalLatitude >= 47.5 & a$decimalLatitude < 48, ]
  cdf48 <- count_species(r48)
  if(nrow(cdf48) > 3){
    shannon_index48<- diversity(cdf48[,-1], index = "shannon")
  }else{
    shannon_index48 <- c("N/A")
  }
  r48.5 <- a[a$decimalLatitude >= 48 & a$decimalLatitude < 48.5, ]
  cdf48.5 <- count_species(r48.5)
  if(nrow(cdf48.5) > 3){
    shannon_index48.5<- diversity(cdf48.5[,-1], index = "shannon")
  }else{
    shannon_index48.5 <- c("N/A")
  }
  
  r49 <- a[a$decimalLatitude >= 48.5 & a$decimalLatitude < 49, ]
  cdf49 <- count_species(r49)
  if(nrow(cdf49) > 3){
    shannon_index49<- diversity(cdf49[,-1], index = "shannon")
  }else{
    shannon_index49 <- c("N/A")
  }
  r49.5 <- a[a$decimalLatitude >= 49 & a$decimalLatitude < 49.5, ]
  cdf49.5 <- count_species(r49.5)
  if(nrow(cdf49.5) > 3){
    shannon_index49.5<- diversity(cdf49.5[,-1], index = "shannon")
  }else{
    shannon_index49.5 <- c("N/A")
  }
  r50 <- a[a$decimalLatitude >= 49.5 & a$decimalLatitude < 50, ]
  cdf50 <- count_species(r50)
  if(nrow(cdf50) > 3){
    shannon_index50<- diversity(cdf50[,-1], index = "shannon")
  }else{
    shannon_index50 <- c("N/A")
  }
  shannon <- na.omit(rbind(shannon_index50, shannon_index49.5, shannon_index49, shannon_index48.5, shannon_index48, shannon_index47.5, shannon_index47, shannon_index46.5, shannon_index46, shannon_index45.5, shannon_index45, shannon_index44.5, shannon_index44, shannon_index43.5, shannon_index43, shannon_index42.5, shannon_index42, shannon_index41.5,  shannon_index41, shannon_index40.5, shannon_index40, shannon_index39.5, shannon_index39, shannon_index38.5, shannon_index38, shannon_index37.5, shannon_index37, shannon_index36.5,  shannon_index36, shannon_index35.5, shannon_index35, shannon_index34.5, shannon_index34, shannon_index33.5, shannon_index33, shannon_index32.5, shannon_index32, shannon_index31.5, shannon_index31, shannon_index30.5))
  col_names <- colnames(shannon)
  col_names[1] <- "Shannon_Diversity_Index"
  colnames(shannon) <- col_names
  shannon <- as.data.frame(shannon)
  shannon <- mutate(data.frame(shannon), lat = c("50", "49.5", "49", "48.5", "48", "47.5", "47", "46.5", "46", "45.5", "45", "44.5", "44", "43.5", "43", "42.5", "42", "41.5",  "41", "40.5", "40", "39.5", "39", "38.5", "38", "37.5", "37", "36.5",  "36", "35.5", "35", "34.5", "34", "33.5", "33", "32.5", "32", "31.5", "31", "30.5"))    
  return(print(shannon))
}

  
df32_36 <- rbind(df32, df34, df36)
df34_36 <- rbind(df34, df36)
df36_38 <- rbind(df36, df38)
df38_40 <- rbind(df38, df40)
df32_34 <- rbind(df32, df34)
df36_38 <- rbind(df36, df38)

# Subsetting the Raw Data
## Mississippi River Forest

LA_Mississippi_River_Forest3 <- less_mlatlongdf(df32, 31.181315, -92.227372, 31.958367, -91.869851, "Louisiana")
LA_Mississippi_River_Forest7 <- ngreater_mlatlongdf(df32_34, 32.470733, -92.293039, 31.958367, -91.869851, "Louisiana")
LA_Mississippi_River_Forest4 <- less_mlatlongdf(df34, 32.470733, -92.293039, 33.053617, -91.694739, "Louisiana")
LA_Mississippi_River_Forest5 <- rbind (LA_Mississippi_River_Forest4, LA_Mississippi_River_Forest3, LA_Mississippi_River_Forest7)
LA_Mississippi_River_Forest6 <- LA_Mississippi_River_Forest5[LA_Mississippi_River_Forest5$decimalLatitude > 31.181315, ]
LA_Mississippi_River_Forest8 <- nless_mlatlongdf(df32_34, 32.746492, -93.748597, 31.181315, -92.227372, "Louisiana")
LA_Mississippi_River_Forest9 <- ngreater_mlatlongdf(LA_Mississippi_River_Forest8, 32.746492, -93.927354, 31.181315, -92.488578, "Louisiana")
LA_Mississippi_River_Forest10 <- ngreater_mlatlongdf(df32, 31.181315, -92.488578, 30.375990, -92.262080, "Louisiana")
LA_Mississippi_River_Forest11 <- ngreater_mlatlongdf(df32, 30.375990, -92.262080, 30.165326, -91.755637, "Louisiana")
LA_Mississippi_River_Forest11 <- LA_Mississippi_River_Forest11[LA_Mississippi_River_Forest11$decimalLongitude < -91.269856, ]
LA_Mississippi_River_Forest12 <- rbind(LA_Mississippi_River_Forest10, LA_Mississippi_River_Forest11)
LA_Mississippi_River_Forest13 <- nless_mlatlongdf(LA_Mississippi_River_Forest12, 31.202187, -91.737765, 30.379166, -90.991478, "Louisiana")
LA_Mississippi_River_Forest14 <- rbind(LA_Mississippi_River_Forest3, LA_Mississippi_River_Forest5, LA_Mississippi_River_Forest9, LA_Mississippi_River_Forest13)
LA_Mississippi_River_Forest15 <- LatLongdf(df34, 33.024550, -93.880469, 32.753200, -93.777361)
LA_Mississippi_River_Forest15 <- LA_Mississippi_River_Forest15[LA_Mississippi_River_Forest15$stateProvince == "Louisiana", ]
LA_Mississippi_River_Forest16 <- rbind(LA_Mississippi_River_Forest14, LA_Mississippi_River_Forest15)
LA_Mississippi_River_Forest <- LA_Mississippi_River_Forest16[!duplicated(LA_Mississippi_River_Forest16$gbifID), ]

LA_MRF_shannon <- gen_shannon(LA_Mississippi_River_Forest)

MI_Mississippi_River_Forest1 <- greater_mlatlongdf(df32_34, 31.998978, -91.109433, 33.687743, -90.024770, "Mississippi")
MI_Mississippi_River_Forest2 <- nless_mlatlongdf(df32_36, 34.428043, -90.096257, 33.687743, -90.024770, "Mississippi")
MI_Mississippi_River_Forest3 <- nless_mlatlongdf(df34_36, 34.541684, -90.229018, 34.428043, -90.096257, "Mississippi")
MI_Mississippi_River_Forest3 <- MI_Mississippi_River_Forest3[MI_Mississippi_River_Forest3$decimalLatitude > 34.428043, ]
MI_Mississippi_River_Forest4 <- LatLongdf(df36, 35.022348, -90.525091, 34.541684, -90.229018)
MI_Mississippi_River_Forest4 <- MI_Mississippi_River_Forest4[MI_Mississippi_River_Forest4$stateProvince == "Mississippi", ]
MI_Mississippi_River_Forest5 <- rbind(MI_Mississippi_River_Forest3, MI_Mississippi_River_Forest1, MI_Mississippi_River_Forest2, MI_Mississippi_River_Forest4)
MI_Mississippi_River_Forest5 <- MI_Mississippi_River_Forest5[!duplicated(MI_Mississippi_River_Forest5$gbifID), ]

MI_MRF_shannon <- gen_shannon(MI_Mississippi_River_Forest5)

AK_Mississippi_River_Forest1 <- less_mlatlongdf(df34, 33.053178, -91.724667, 33.605452, -91.534211, "Arkansas")
AK_Mississippi_River_Forest2 <- ngreater_mlatlongdf(df34_36, 34.747772, -92.176045, 33.605452, -91.534211, "Arkansas")
AK_Mississippi_River_Forest2 <- AK_Mississippi_River_Forest2[AK_Mississippi_River_Forest2$decimalLatitude > 33.605452, ]
AK_Mississippi_River_Forest3 <- less_mlatlongdf(df36_38, 34.747772, -92.176045, 36.526765, -90.732280, "Arkansas")
AK_Mississippi_River_Forest3 <- AK_Mississippi_River_Forest3[AK_Mississippi_River_Forest3$decimalLatitude > 34.747772, ]
AK_Mississippi_River_Forest4 <- LatLongdf(df34, 33.438286, -93.797359, 33.050747, -93.747094)
AK_Mississippi_River_Forest <- rbind(AK_Mississippi_River_Forest1, AK_Mississippi_River_Forest2, AK_Mississippi_River_Forest3, AK_Mississippi_River_Forest4)
AK_Mississippi_River_Forest <- AK_Mississippi_River_Forest[!duplicated(AK_Mississippi_River_Forest$gbifID), ]

AK_MRF_shannon <- gen_shannon(AK_Mississippi_River_Forest)

MO_Mississippi_River_Forest <- less_mlatlongdf(df38, 36.509215, -90.736273, 37.285163, -89.787756, "Missouri")

MO_MRF_shannon <- gen_shannon(MO_Mississippi_River_Forest)

TN_Mississippi_River_Forest <- greater_mlatlongdf(df36_38, 34.997299, -90.139402, 36.510921, -89.330913, "Tennessee")

TN_MRF_shannon <- gen_shannon(TN_Mississippi_River_Forest)

KY_Mississippi_River_Forest <- greater_mlatlongdf(df38, 36.493387, -89.219074, 36.950273, -89.101720, "Kentucky")

KY_MRF_Shannon <- gen_shannon(KY_Mississippi_River_Forest)

Mississippi_River_Forest <- rbind(AK_Mississippi_River_Forest, KY_Mississippi_River_Forest, TN_Mississippi_River_Forest, LA_Mississippi_River_Forest, MO_Mississippi_River_Forest, MI_Mississippi_River_Forest5)

## Mixed Forests

MI_Mixed_Forest1 <- greater_mlatlongdf(df32_34, 31.029375, -90.670486, 32.035643, -89.635176, "Mississippi")
MI_Mixed_Forest2 <- greater_mlatlongdf(df34_36, 32.035643, -89.635176, 32.542343, -88.454560, "Mississippi")
MI_Mixed_Forest3 <- rbind(MI_Mixed_Forest1, MI_Mixed_Forest2)
MI_Mixed_Forest3 <- MI_Mixed_Forest3[!duplicated(MI_Mixed_Forest3$gbifID), ]
MI_Mixed_Forest4 <- anti_join(MI_Mixed_Forest3, MI_Mississippi_River_Forest5)

MI_MF_shannon <- gen_shannon(MI_Mixed_Forest4)

SC_Mixed_Forest1<- greater_mlatlongdf(df34_36, 33.288275, -81.863640, 34.819184, -80.575592, "South Carolina")
SC_Mixed_Forest <- less_mlatlongdf(SC_Mixed_Forest1, 33.288275, -82.899328, 35.197685, -82.226109, "South Carolina")

SC_MF_shannon <- gen_shannon(SC_Mixed_Forest)

VA_Mixed_Forest1 <- greater_mlatlongdf(df38_40, 36.572383, -77.667126, 38.502839, -77.304794, "Virginia")
VA_Mixed_Forest <- less_mlatlongdf(VA_Mixed_Forest1, 36.560859, -80.987262, 39.217619, -77.455949, "Virginia")

VA_MF_shannon <- gen_shannon(VA_Mixed_Forest)

NC_Mixed_Forest1 <- greater_mlatlongdf(df36_38, 34.796832, -79.636560, 36.559827, -77.722891, "North Carolina")
NC_Mixed_Forest <- less_mlatlongdf(NC_Mixed_Forest1, 34.796832, -82.370233, 36.562920, -80.741205, "North Carolina")

NC_MF_shannon <- gen_shannon(NC_Mixed_Forest)

GA_Mixed_Forest1 <- greater_mlatlongdf(df32_36, 31.406694, -85.068045, 33.288275, -81.863640, "Georgia")
GA_Mixed_Forest <- GA_Mixed_Forest1[GA_Mixed_Forest1$decimalLatitude < 34.456300, ]

GA_MF_shannon <- gen_shannon(GA_Mixed_Forest)

TX_Mixed_Forest3 <- nless_mlatlongdf(df32_34, 30.488458, -94.396987, 30.000992, -93.799283, "Texas")
TX_Mixed_Forest4 <- TX_Mixed_Forest3[TX_Mixed_Forest3$decimalLongitude > -94.353667, ]
TX_Mixed_Forest5 <- greater_mlatlongdf(df32_34, 31.833599, -93.884635, 30.810800, -94.438310, "Texas")
TX_Mixed_Forest6 <- TX_Mixed_Forest5[TX_Mixed_Forest5$decimalLongitude > -94.438310, ]
TX_Mixed_Forest7 <- rbind(TX_Mixed_Forest4, TX_Mixed_Forest6)
TX_Mixed_Forest9 <- TX_Mixed_Forest7[!duplicated(TX_Mixed_Forest7$gbifID), ]

TX_MF_shannon <- gen_shannon(TX_Mixed_Forest9)

LA_Mixed_Forest1 <- greater_mlatlongdf(df32_34, 31.711728, -93.830905, 32.373758, -92.178042, "Louisiana")
LA_Mixed_Forest <- anti_join(LA_Mixed_Forest1, LA_Mississippi_River_Forest)

LA_MF_shannon <- gen_shannon(LA_Mixed_Forest)

AK_Mixed_Forest1 <- LatLongdf(df34_36, 34.081595, -94.457706, 33.578257, -93.176890)
AK_Mixed_Forest2 <- less_mlatlongdf(df34_36, 34.046414, -93.202459, 34.690341, -92.281856, "Arkansas")
AK_Mixed_Forest3 <- AK_Mixed_Forest2[AK_Mixed_Forest2$decimalLongitude < -93.202459, ]
AK_Mixed_Forest4 <- LatLongdf(df34_36, 33.910481, -93.202459, 33.259865, -91.697607)
AK_Mixed_Forest5 <- rbind(AK_Mixed_Forest3, AK_Mixed_Forest4)
AK_Mississippi_River_Forest1 <- LatLongdf(df34_36, 33.537548, -93.902493, 33.028345, -93.656557)
AK_Mixed_Forest7 <- anti_join(AK_Mixed_Forest5, AK_Mississippi_River_Forest1)
AK_Mixed_Forest8 <- LatLongdf(df34_36, 35.580192, -94.404706, 35.081387, -92.272595)
AK_Mixed_Forest <- rbind(AK_Mixed_Forest7, AK_Mixed_Forest8)

AK_MF_shannon <- gen_shannon(AK_Mixed_Forest)

AL_Mixed_Forest1 <- LatLongdf(df36, 34.456300, -88.157066, 31.547402, -84.588848)
AL_Mixed_Forest2 <- LatLongdf(df34, 34.456300, -88.157066, 31.547402, -84.588848)
AL_Mixed_Forest3 <- LatLongdf(df32, 34.456300, -88.157066, 31.547402, -84.588848)
AL_Mixed_Forest4 <- rbind(AL_Mixed_Forest1, AL_Mixed_Forest2, AL_Mixed_Forest3)
AL_Mixed_Forest <- AL_Mixed_Forest4[AL_Mixed_Forest4$stateProvince == "Alabama" & AL_Mixed_Forest4$decimalLatitude >= 31.433183, ]

AL_MF_shannon <- gen_shannon(AL_Mixed_Forest)

Mixed_Forest <- rbind(TX_Mixed_Forest9, AL_Mixed_Forest, LA_Mixed_Forest, AK_Mixed_Forest, GA_Mixed_Forest, MI_Mixed_Forest4, SC_Mixed_Forest, NC_Mixed_Forest, VA_Mixed_Forest)
## Coastal Mixed Forests

Louisiana <- df32_34[df32_34$stateProvince == "Louisiana", ]
LA_Coastal_Mixed_Forest1 <- anti_join(Louisiana, LA_Mixed_Forest)
LA_Coastal_Mixed_Forest2 <- anti_join(LA_Coastal_Mixed_Forest1, LA_Mississippi_River_Forest)
LA_Coastal_Mixed_Forest3 <- LA_Coastal_Mixed_Forest2[LA_Coastal_Mixed_Forest2$decimalLatitude < 32.081824, ]
LA_Coastal_Mixed_Forest4 <- draw_complex(30.671334, -91.358859,31.036490, -91.015832, "greater", 31.036490, -91.015832, 31.012974, -91.633281, "less", 30.671334, -91.358859, 31.012974, -91.633281, "ngreater", "Louisiana")
LA_Coastal_Mixed_Forest5 <- undup(LA_Coastal_Mixed_Forest4)
LA_Coastal_Mixed_Forest <- anti_join(LA_Coastal_Mixed_Forest3, LA_Coastal_Mixed_Forest5)

LA_CMF_shannon <- gen_shannon(LA_Coastal_Mixed_Forest)

MI_Coastal_Mixed_Forest1 <- less_mlatlongdf(df32_34, 31.012218, -89.794932, 31.540894, -89.420732, "Mississippi")
MI_Coastal_Mixed_Forest2 <- MI_Coastal_Mixed_Forest1[MI_Coastal_Mixed_Forest1$decimalLatitude < 31.547402, ]
MI_Coastal_Mixed_Forest3 <- lat_less(df32, 31.012218, "Mississippi")
MI_Coastal_Mixed_Forest <- rbind(MI_Coastal_Mixed_Forest2, MI_Coastal_Mixed_Forest3)

MI_CMF_shannon <- gen_shannon(MI_Coastal_Mixed_Forest)

Alabama <- df32_36[df32_36$decimalLatitude < 31.547402 & df32_36$stateProvince == "Alabama", ]
AL_Coastal_Mixed_Forest <- anti_join(Alabama, AL_Mixed_Forest)

AL_CMF_shannon <- gen_shannon(AL_Coastal_Mixed_Forest)

NC_Coastal_Mixed_Forest1 <- draw_complex(35.190545, -78.548216, 34.647278, -79.523298, "less", 35.905851, -77.681936, 35.190545, -78.548216, "less", 36.592608, -77.403379, 35.905851, -77.681936, "less", "North Carolina")
NC_Coastal_Mixed_Forest2 <- lat_less(df34_36, 34.647278, "North Carolina")
NC_Coastal_Mixed_Forest <- rbind (NC_Coastal_Mixed_Forest1, NC_Coastal_Mixed_Forest2)

NC_CMF_shannon <- gen_shannon(NC_Coastal_Mixed_Forest)

SC_Coastal_Mixed_Forest1 <- less_mlatlongdf(df34_36, 33.414474, -82.011009, 34.832187, -79.560801, "South Carolina")
SC_Coastal_Mixed_Forest2 <- lat_less(df34, 33.414474, "South Carolina")
SC_Coastal_Mixed_Forest <- rbind(SC_Coastal_Mixed_Forest1, SC_Coastal_Mixed_Forest2)

SC_CMF_shannon <- gen_shannon(SC_Coastal_Mixed_Forest)

VA_Coastal_Mixed_Forest <- draw_complex(37.164666, -77.334652, 36.521970, -77.575582, "less", 37.164666, -77.334652, 37.907627, -77.401577, "ngreater", 39.046219, -77.243896, 37.907627, -77.401577, "less", "Virginia")

VA_CMF_shannon <- gen_shannon(VA_Coastal_Mixed_Forest)

GA_Coastal_Mixed_Forest1 <- less_mlatlongdf(df32_34, 31.685326, -85.206560, 33.329061, -81.826020, "Georgia")
GA_Coastal_Mixed_Forest2 <- lat_less(df32, 31.685326, "Georgia")
GA_Coastal_Mixed_Forest <- rbind(GA_Coastal_Mixed_Forest1, GA_Coastal_Mixed_Forest2)

GA_CMF_shannon <- gen_shannon(GA_Coastal_Mixed_Forest)

MD_Coastal_Mixed_Forest1 <- draw_complex(39.370324, -76.859077, 39.038422, -77.304129, "less", 39.513820, -76.238534, 39.370324, -76.859077, "less", 39.513820, -76.238534, 39.513820, -76.238534, "less", "Maryland" )                         
MD_Coastal_Mixed_Forest2 <- lat_less(df40, 39.038422, "Maryland")
MD_Coastal_Mixed_Forest <- rbind(MD_Coastal_Mixed_Forest1, MD_Coastal_Mixed_Forest2)

MD_CMF_shannon <- gen_shannon(MD_Coastal_Mixed_Forest)

DA_Coastal_Mixed_Forest1 <- lat_less(df40, 39.158189, "Delaware")
DA_Coastal_Mixed_Forest2 <- nless_mlatlongdf(df40, 39.492950, -75.790751, 39.158189, -75.382036, "Delaware")
DA_Coastal_Mixed_Forest <- rbind(DA_Coastal_Mixed_Forest1, DA_Coastal_Mixed_Forest2)

DA_CMF_shannon <- gen_shannon(DA_Coastal_Mixed_Forest)

Coastal_Mixed_Forest <- rbind(MD_Coastal_Mixed_Forest, DA_Coastal_Mixed_Forest, MI_Coastal_Mixed_Forest, GA_Coastal_Mixed_Forest, AL_Coastal_Mixed_Forest, MI_Coastal_Mixed_Forest, LA_Coastal_Mixed_Forest, SC_Coastal_Mixed_Forest, NC_Coastal_Mixed_Forest, VA_Coastal_Mixed_Forest)

## Appalachian Broadleaf Forest

GA_ABF1 <- greater_mlatlongdf(df36, 34.516167, -84.521330, 34.999213, -83.195299, "Georgia")
GA_ABF2 <- ngreater_mlatlongdf(GA_ABF1, 35.011333, -84.692971, 34.516167, -84.521330, "Georgia")

GA_ABF_shannon <- gen_shannon(GA_ABF2)

TN_ABF1 <- less_mlatlongdf(df36, 35.014611, -84.658249, 35.449475, -84.379607, "Tennessee")
TN_ABF2 <- less_mlatlongdf(df36, 35.449475, -84.379607, 35.864233, -83.205525, "Tennessee")
TN_ABF3 <- less_mlatlongdf(df36_38, 35.864233, -83.205525, 36.206725, -82.451103, "Tennessee")
TN_ABF4 <- less_mlatlongdf(df38, 36.206725, -82.451103, 36.457293, -82.166540, "Tennessee")
TN_ABF5 <- less_mlatlongdf(df38, 36.457293, -82.166540, 36.632741, -81.762858, "Tennessee")
TN_ABF <- rbind(TN_ABF1, TN_ABF2, TN_ABF3, TN_ABF4, TN_ABF5)
TN_ABF <- TN_ABF[!duplicated(TN_ABF$gbifID), ]

TN_ABF_shannon <- gen_shannon(TN_ABF)

NC_ABF1 <- greater_mlatlongdf(df36, 35.174796, -82.576840, 35.654798, -82.219482, "North Carolina")
NC_ABF2 <- greater_mlatlongdf(df36_38, 35.654798, -82.219482, 36.594499, -80.877589, "North Carolina")
NC_ABF <- rbind(NC_ABF1, NC_ABF2)

NC_ABF_shannon <- gen_shannon(NC_ABF)

df38_40 <- rbind(df38, df40)
VA_ABF1 <- greater_mlatlongdf(df38, 36.594499, -80.877589, 36.971052, -80.615299, "Virginia")
VA_ABF2 <- greater_mlatlongdf(df38, 36.971052, -80.615299, 37.869152, -79.156469, "Virginia")
VA_ABF3 <- greater_mlatlongdf(df38_40, 37.869152, -79.156469, 39.213878, -77.930084, "Virginia")                             
VA_ABF <- rbind(VA_ABF1, VA_ABF2, VA_ABF3)
VA_ABF <- VA_ABF[VA_ABF$decimalLongitude > -82.513925, ]

VA_ABF_shannon <- gen_shannon(VA_ABF)


WV_ABF1 <- less_mlatlongdf(df38_40, 37.539893, -82.000615, 39.749058, -79.787681, "West Virginia")
WV_ABF2 <- ngreater_mlatlongdf(df38, 37.539893, -82.000615, 37.156101, -81.579996, "West Virginia")
WV_ABF <- rbind(WV_ABF1, WV_ABF2)


WV_ABF_shannon <- gen_shannon(WV_ABF)

df40_42 <- rbind(df40, df42)
MD_ABF <- greater_mlatlongdf(df40_42, 39.386051, -77.777937, 39.799047, -77.446161, "Maryland")

MD_ABF_shannon <- gen_shannon(MD_ABF)

PA_ABF1 <- greater_mlatlongdf(df40_42, 39.699045, -77.523769, 40.152996, -77.156800, "Pennsylvania")
PA_ABF2 <- greater_mlatlongdf(df40_42, 40.152996, -77.156800, 40.577551, -76.123633, "Pennsylvania")
PA_ABF3 <- nless_mlatlongdf(df42, 41.148202, -76.518326, 40.577551, -76.123633, "Pennsylvania")
PA_ABF4 <- less_mlatlongdf(df40_42, 41.078234, -77.505058, 41.148202, -76.518326, "Pennsylvania")
PA_ABF5 <- less_mlatlongdf(df40_42, 39.707910, -79.652651, 41.078234, -77.505058, "Pennsylvania")
PA_ABF6 <- rbind(PA_ABF1, PA_ABF2, PA_ABF3)
PA_ABF7 <- rbind(PA_ABF4, PA_ABF5)
PA_ABF <- intersect(PA_ABF6, PA_ABF7)

PA_ABF_shannon <- gen_shannon(PA_ABF)

SC_ABF <- greater_mlatlongdf(df36, 34.605987, -83.232190, 35.251023, -82.312916, "South Carolina")

SC_ABF_shannon <- gen_shannon(SC_ABF)

Appalacian_Broadleaf_Forest <- rbind(GA_ABF2, SC_ABF, NC_ABF, VA_ABF, MD_ABF, PA_ABF, WV_ABF, TN_ABF)

## Adirondack
df42_44 <- rbind(df42, df44)
NY_A1 <- greater_mlatlongdf(df42_44, 42.337786, -75.689335, 42.374406, -74.731311, "New York")
NY_A2 <- greater_mlatlongdf(df42_44, 42.374406, -74.731311, 42.669318, -74.702730, "New York")
NY_A3 <- nless_mlatlongdf(df42_44, 42.852931, -75.288654, 42.669318, -74.702730, "New York")
NY_A4 <- less_mlatlongdf(df42_44, 42.464080, -75.531598, 42.852931, -75.288654, "New York")
NY_A5 <- less_mlatlongdf(df42_44, 42.337786, -75.689335, 42.464080, -75.531598, "New York")
NY_A6 <- rbind(NY_A1, NY_A2)
NY_A6 <- NY_A6[NY_A6$decimalLongitude > -75.689335, ]
NY_A7 <- rbind(NY_A3, NY_A4, NY_A5)
NY_A7 <- NY_A7[NY_A7$decimalLongitude < -74.702730, ]
NY_A8 <- rbind(NY_A7, NY_A6)
NY_A8 <- NY_A8[!duplicated(NY_A8$gbifID), ]


df44_46 <- rbind(df44, df46)
NY_A9 <- LatLongdf(df44_46, 44.559450, -74.883355, 43.496556, -74.009059)
NY_A <- rbind(NY_A8, NY_A9)
NY_A_shannon <- gen_shannon(NY_A)

df46_48 <- rbind(df46, df48)
NE_A1 <- LatLongdf(df44_46, 44.997276, -72.730665, 42.734387, -71.726884)
NE_A2 <- greater_mlatlongdf(df44, 42.734387, -71.726884, 43.803408, -71.634483, "New Hampshire")
NE_A3 <- greater_mlatlongdf(df44_46, 43.803408, -71.634483, 44.211874, -71.031624, "New Hampshire")
NE_A4 <- greater_mlatlongdf(df46, 44.211874, -71.031624, 45.489656, -69.499496, "Maine")
NE_A5 <- greater_mlatlongdf(df46_48, 45.489656, -69.499496, 47.507471, -69.086041, "Maine")
NE_A6 <- LatLongdf(df44, 42.760428, -73.435186, 42.041722, -72.825415) 
NE_A7 <- df44_46[df44_46$decimalLatitude > 44.211874 & df44_46$stateProvince == "New Hampshire", ]
NE_A <- rbind(NE_A1, NE_A2, NE_A3, NE_A4, NE_A5, NE_A6, NE_A7)
NE_A <- NE_A[!duplicated(NE_A$gbifID), ]
NE_A_shannon <- gen_shannon(NE_A)

Adirondack <- rbind(NE_A, NY_A)

## Laurentian Mixed Forest

PA_LMF1 <- less_mlatlongdf(df42, 41.476261, -79.674768, 42.051306, -79.050664, "Pennsylvania")
PA_LMF2 <- ngreater_mlatlongdf(df42, 41.476261, -79.674768, 41.134510, -77.819871, "Pennsylvania")
PA_LMF2 <- PA_LMF2[PA_LMF2$decimalLongitude > -79.674768, ]
PA_LMF3 <- ngreater_mlatlongdf(df42, 41.476261, -79.674768, 41.462520, -76.098071, "Pennsylvania")
PA_LMF4 <- greater_mlatlongdf(df42, 41.476261, -79.674768, 41.760609, -75.505725, "Pennsylvania")
PA_LMF4 <- PA_LMF4[PA_LMF4$decimalLongitude > -79.674768, ]
PA_LMF5 <- ngreater_mlatlongdf(df42, 41.760609, -75.505725, 41.291566, -75.266508, "Pennsylvania")
PA_LMF6 <- ngreater_mlatlongdf(df42, 41.291566, -75.266508, 41.036971, -75.011659, "Pennsylvania")
PA_LMF7 <- rbind (PA_LMF1, PA_LMF2, PA_LMF3, PA_LMF4)
PA_LMF8 <- intersect(PA_LMF7, PA_ABF)
PA_LMF9 <- anti_join(PA_LMF7, PA_LMF8)
PA_LMF <- rbind(PA_LMF9, PA_LMF5, PA_LMF6)
PA_LMF_shannon <- gen_shannon(PA_LMF)

NY_LMF1 <- LatLongdf(df44, 42.790683, -78.823272, 42.079885, -76.439192)
NY_LMF2 <- less_mlatlongdf(df44, 42.790683, -78.823272, 43.033518, -75.097159, "New York")
NY_LMF3 <- rbind(NY_LMF1, NY_LMF2)
NY_LMF4 <- anti_join(NY_LMF3, NY_A)
NY_LMF <- NY_LMF4[NY_LMF4$decimalLongitude < -74.181420, ]

NY_LMF_shannon <- gen_shannon(NY_LMF)

df46_48 <- rbind(df46, df48)
df44_46 <- rbind(df44, df46)
ME_LMF1 <- LatLongdf(df46_48, 47.315615, -68.481730, 44.346342, -66.934446)
ME_LMF2 <- less_mlatlongdf(df44_46, 43.607842, -71.023697, 45.298208, -68.323844, "Maine")
ME_LMF2 <- ME_LMF2[ME_LMF2$decimalLongitude > -69.762756 & ME_LMF2$stateProvince == "Maine", ]
ME_LMF <- rbind (ME_LMF1, ME_LMF2)
ME_LMF <- ME_LMF[!duplicated(ME_LMF$gbifID), ]

ME_LMF_shannon <- gen_shannon(ME_LMF)

MN_LMF <- ngreater_mlatlongdf(df48, 48.692677, -94.498492, 46.059647, -92.318089, "Minnesota")

MN_LMF_shannon <- gen_shannon(MN_LMF)

WI_LMF1 <- ngreater_mlatlongdf(df46_48, 46.060251, -92.439432,45.300294, -87.640502, "Wisconsin")
WI_LMF2 <- df46_48[df46_48$decimalLatitude > 46.060251 & df46_48$stateProvince == "Wisconsin", ]
WI_LMF <- rbind(WI_LMF1, WI_LMF2)

WI_LMF_shannon <- gen_shannon(WI_LMF)

MI_LMF1 <- greater_mlatlongdf(df46_48, 44.301494, -86.382115, 44.833314, -82.777586, "Michigan")
MI_LMF2 <- df46_48[df46_48$decimalLatitude > 44.833314 & df46_48$stateProvince == "Michigan", ]
MI_LMF <- rbind(MI_LMF1, MI_LMF2)

MI_LMF_shannon <- gen_shannon(MI_LMF)

LMF <- rbind(PA_LMF, NY_LMF, ME_LMF, MI_LMF, MN_LMF, WI_LMF)

## Eastern Broadleaf Forest
df42_44<-rbind(df42,df44)
PA_EBF1 <- LatLongdf(df42_44, 42.219068, -80.524737, 40.106398, -79.901035)
PA_EBF2 <- greater_mlatlongdf(df42, 40.106398, -79.901035, 40.936084, -79.099133, "Pennsylvania")
PA_EBF3 <- nless_mlatlongdf(df42, 41.296883, -79.867623, 40.106398, -79.901035, "Pennsylvania")
PA_EBF4 <- rbind(PA_EBF1, PA_EBF2, PA_EBF3)
PA_EBF4 <- PA_EBF4[!duplicated(PA_EBF4$gbifID), ]

PA_EBF5 <- less_mlatlongdf(df40_42, 39.687710, -76.559776, 40.505590, -75.033934, "Pennsylvania")
PA_EBF <- rbind(PA_EBF4, PA_EBF5)

PA_EBF_shannon <- gen_shannon(PA_EBF)

NE_EBF1 <- LatLongdf(df42, 41.708336, -73.751449, 41.342060, -71.825478)
RI <- df42_44[df42_44$stateProvince == "Rhode Island", ]
NJ <- df40_42[df40_42$stateProvince == "New Jersey", ]
MA_EBF <- df44[df44$stateProvince == "Massachusetts" & df44$decimalLongitude > -71.685720, ]

MA_EBF_shannon <- gen_shannon(MA_EBF)

TN_EBF1 <- less_mlatlongdf(df36_38, 34.965046, -86.239315, 36.717117, -84.787055, "Tennessee")
TN_EBF2 <- greater_mlatlongdf(TN_EBF1, 34.923520, -84.922149, 35.461722, -84.668848, "Tennessee")
TN_EBF3 <- greater_mlatlongdf(TN_EBF1, 35.461722, -84.668848, 36.635856, -82.726873, "Tennessee")
TN_EBF <- rbind(TN_EBF2, TN_EBF3)

TN_EBF_shannon <- gen_shannon(TN_EBF)

KY_EBF1<- less_mlatlongdf(df38_40, 36.595193, -85.040356, 38.760076, -83.351682, "Kentucky")
KY_EBF<- greater_mlatlongdf(KY_EBF1, 36.595193, -84.111585, 38.937563, -81.859382, "Kentucky")

KY_EBF_shannon <- gen_shannon(KY_EBF)

WV_EBF <- greater_mlatlongdf(df40_42, 38.018850, -82.566057, 40.624830, -79.247208, "West Virginia")

WV_EBF_shannon <- gen_shannon(WV_EBF)

OH_EBF <- less_mlatlongdf(df40_42, 38.306586, -83.790119, 40.854302, -80.446032, "Ohio")

OH_EBF_shannon <- gen_shannon(OH_EBF)

NH_EBF <- less_mlatlongdf(df44, 42.736690, -71.521504, 43.345396, -70.890544, "New Hampshire")

NH_EBF_shannon <- gen_shannon(NH_EBF)

ME_EBF <- less_mlatlongdf(df44, 43.345396, -70.890544, 43.437097, -70.347918, "Maine")

ME_EBF_shannon <- gen_shannon(ME_EBF)

EBF <- rbind(ME_EBF, NH_EBF, OH_EBF, WV_EBF, KY_EBF, TN_EBF, NE_EBF1, MA_EBF, RI, NJ, PA_EBF)

## Eastern Broadleaf Forest (Continental)

TN_EBFC1 <- greater_mlatlongdf(df36_38, 34.941742, -87.262651, 36.738155, -85.594801, "Tennessee")
TN_EBFC2 <- LatLongdf(df36, 35.825465, -89.910586, 34.987008, -87.188437)
TN_EBFC3 <- greater_mlatlongdf(TN_EBFC1, 34.999973, -90.116330, 36.578941, -89.214222, "Tennessee")
TN_EBFC4 <- anti_join(TN_EBFC1, TN_EBFC2)
TN_EBFC5 <- anti_join(TN_EBFC4, TN_EBFC3)
TN_EBFC6 <- nless_mlatlongdf(df36, 35.594151, -90.037198,34.999973, -89.483272, "Tennessee")
TN_EBFC <- anti_join(TN_EBFC5, TN_EBFC6)

KY_EBFC1 <- greater_mlatlongdf(df38_40, 36.518129, -85.823274, 38.939155, -83.926950, "Kentucky")
KY_EBFC2 <- long_less(df38, -89.008478, "Kentucky")
KY_EBFC <- anti_join(KY_EBFC1, KY_EBFC2)

OH_EBFC <- greater_mlatlongdf(df40_42, 38.689926, -84.018339, 41.605665, -82.099169, "Ohio")

IN_EBFC1 <- less_mlatlongdf(df40_42, 39.292184, -87.626524, 41.878436, -84.681071, "Indiana")
IN_EBFC2 <- df40_42[df40_42$decimalLatitude < 39.292184 & df40_42$stateProvince == "Indiana", ]
IN_EBFC <- rbind(IN_EBFC1, IN_EBFC2)

IL_EBFC1 <- df38_40[df38_40$decimalLatitude < 38.795748 & df38_40$decimalLatitude > 37.731468 & df38_40$stateProvince == "Illinois", ]
IL_EBFC2 <- ngreater_mlatlongdf(df42_44, 42.541959, -88.670219, 41.769511, -87.562259, "Illinois")              
IL_EBFC <- rbind(IL_EBFC1, IL_EBFC2)

MO_EBFC1 <- less_mlatlongdf(df38_40, 36.901188, -94.723119, 38.455479, -90.263935, "Missouri")
MO_EBFC3 <- df38_40[df38_40$decimalLatitude<36.901188 & df38_40$decimalLatitude > 36.498232 & df38_40$stateProvince == "Missouri", ]
MO_EBFC4 <- rbind(MO_EBFC1, MO_EBFC3)
MO_EBFC2 <- less_mlatlongdf(MO_EBFC4, 36.475191, -91.052744, 37.733348, -89.475126, "Missouri")
MO_EBFC <- anti_join(MO_EBFC4, MO_EBFC2)

IA_EBFC1 <- lat_greater(df44, 42.356188, "Iowa")
IA_EBFC <- ngreater_mlatlongdf(IA_EBFC1, 43.590333, -93.499367, 42.522129, -92.894690, "Iowa")

MN_EBFC1 <- LatLongdf(df44_46, 44.740897, -92.841155, 43.545701, -91.175522)
MN_EBFC <- MN_EBFC1[MN_EBFC1$stateProvince == "Minnesota", ]

WI_EBFC1 <- nless_mlatlongdf(df44_46, 44.813037, -92.842505, 43.270794, -87.826464, "Wisconsin")
WI_EBFC3 <- df44_46[df44_46$decimalLatitude < 43.270794 & df44_46$stateProvince == "Wisconsin", ]
WI_EBFC4 <- rbind(WI_EBFC1, WI_EBFC3)
WI_EBFC2 <- LatLongdf(df44, 42.938808, -89.702150, 42.468908, -88.934824)
WI_EBFC <- anti_join(WI_EBFC4, WI_EBFC2)

MI_EBFC <- lat_less(df42_44, 42.915728, "Michigan")

AR_EBFC1 <- lat_greater(df38, 36.043312, "Arkansas")
AR_EBFC <- greater_mlatlongdf(AR_EBFC1, 35.838297, -91.651495, 36.642461, -90.970958, "Arkansas")

EBFC <- rbind(AR_EBFC, KY_EBFC, MI_EBFC, TN_EBFC, WI_EBFC, IN_EBFC, IL_EBFC, MI_EBFC, OH_EBFC, MO_EBFC, IA_EBFC, MN_EBFC)

## Temperate Prairie Parkland

MN_TPP <- nless_mlatlongdf(df44_46, 45.998736, -96.755394, 43.590333, -93.499367, "Minnesota")

IA_TPP1 <- nless_mlatlongdf(df42_44, 43.590333, -93.499367, 40.583128, -92.296986, "Iowa")                            
IA_TPP2 <- lat_less(df42_44, 42.241861, "Iowa")
IA_TPP3 <- rbind(IA_TPP1, IA_TPP2)
IA_TPP <- undup(IA_TPP3)

df38_42 <- rbind(df38, df42)
MO_TPP1 <- greater_mlatlongdf(df38_42, 37.577860, -94.766580, 40.721161, -92.305797, "Missouri")
MO_TPP2 <- lat_greater(df40_42, 39.809144, "Missouri")
MO_TPP3 <- rbind(MO_TPP1, MO_TPP2)
MO_TPP <- undup(MO_TPP3)

IL_TPP1 <- greater_mlatlongdf(df40_42, 38.994239, -90.735583, 40.275630, -87.524847, "Illinois")
IL_TPP5 <- lat_greater(df42_44, 40.275630, "Illinois")
IL_TPP6 <- rbind(IL_TPP1, IL_TPP5)
IL_TPP2 <- ngreater_mlatlongdf(df42_44, 42.611734, -89.305985, 41.146083, -87.524847, "Illinois")
IL_TPP3 <- greater_mlatlongdf(df42_44, 41.987715, -90.313735, 42.559969, -89.540346, "Illinois")
IL_TPP4 <- rbind(IL_TPP2, IL_TPP3)
IL_TPP <- anti_join(IL_TPP6, IL_TPP4)

ND_TPP <- ngreater_mlatlongdf(df46_48, 49.057345, -97.733164, 45.908338, -97.424144, "North Dakota")
df44_48<-rbind(df44,df48)
SD_TPP <- ngreater_mlatlongdf(df44_48, 46.028289, -97.377730, 42.765704, -97.206595, "South Dakota")
df40_44 <- rbind(df40, df44)
NE_TPP <- long_greater(df40_44, -97.843503, "Nebraska")

KS_TPP <- long_greater(df38_40, -96.956059, "Kansas")

OK_TPP1 <- greater_mlatlongdf(df38, 36.206196, -95.403032, 36.907900, -94.598785, "Oklahoma")
OK_TPP1 <- OK_TPP1[OK_TPP1$decimalLongitude > -95.403032, ]
OK_TPP2 <- ngreater_mlatlongdf(df38, 37.018695, -96.221144, 36.206196, -95.403032, "Oklahoma")
OK_TPP3 <- LatLongdf(df38, 37.018695,-95.403032, 36.206196, -95.403032)
OK_TPP4 <- OK_TPP3[OK_TPP3$stateProvince == "Oklahoma",]
OK_TPP <- rbind(OK_TPP1, OK_TPP2, OK_TPP4)

TPP <- rbind(OK_TPP, KS_TPP, SD_TPP, ND_TPP, IL_TPP, MO_TPP, IA_TPP, MN_TPP, NE_TPP)

## Great Plains Steppe

ND_GPS1 <- nless_mlatlongdf(df48, 49.047857, -103.653474, 47.027358, -100.340218, "North Dakota")
ND_GPS2 <- nless_mlatlongdf(df46_48, 47.027358, -100.340218, 45.948007, -100.291731, "North Dakota")
ND_GPS3 <- rbind(ND_GPS1, ND_GPS2)
ND_GPS <- long_less(ND_GPS3, -98.125993, "North Dakota")

SD_GPS1 <- greater_mlatlongdf(df44_48, 42.968319, -101.272126, 46.021696, -99.764746, "South Dakota")
SD_GPS2 <- LatLongdf(df44_46, 44.734438, -104.124933, 42.997714, -103.264270)
SD_GPS <- anti_join(SD_GPS1, SD_GPS2)

NE_GPS1 <- long_less(df42_44, -103.079814, "Nebraska")
NE_GPS2 <- nless_mlatlongdf(df40_42, 41.835591, -102.902051, 40.0000, -100.033768, "Nebraska")
NE_GPS3 <- rbind(NE_GPS1, NE_GPS2)
NE_GPS <- undup(NE_GPS3)

KS_GPS <- long_less(df38_40, -100.033768, "Kansas")

OK_GPS <- long_less(df38, -100.078476, "Oklahoma")

CO_GPS1 <- long_greater(df38, -104.020021, "Colorado")
CO_GPS2 <- long_greater(df40_42, -104.020021, "Colorado")
CO_GPS <- rbind(CO_GPS1, CO_GPS2)

WY_GPS1 <- ngreater_mlatlongdf(df42_44, 45.210833, -106.940537, 41.056021, -104.596913, "Wyoming")
WY_GPS2 <- ngreater_mlatlongdf(df46, 45.210833, -106.940537, 41.056021, -104.596913, "Wyoming")
WY_GPS3 <- rbind (WY_GPS1, WY_GPS2)
WY_GPS4 <- LatLongdf(df44_46, 45.026754, -105.026135, 43.642341, -104.007568)
WY_GPS <- anti_join(WY_GPS3, WY_GPS4)

MT_GPS1 <- ngreater_mlatlongdf(df46, 45.782291, -108.443955, 44.026411, -104.045762, "Montana")
MT_GPS2 <- less_mlatlongdf(df46_48, 45.782291, -108.443955, 46.494697, -106.530649, "Montana")
MT_GPS3 <- ngreater_mlatlongdf(df46_48, 49.088501, -111.914258, 46.494697, -106.530649, "Montana")
MT_GPS <- rbind(MT_GPS1, MT_GPS2, MT_GPS3)

GPS <- rbind(ND_GPS, SD_GPS, NE_GPS, KS_GPS, OK_GPS, CO_GPS, WY_GPS, MT_GPS)

## Black Hills Coniferous Forest Pine

SD_BHCF1 <- nless_mlatlongdf(df46, 44.334871, -103.445674, 44.129708, -103.303781, "South Dakota")
SD_BHCF2 <- nless_mlatlongdf(df46,  44.445140, -103.624872, 44.334871, -103.445674, "South Dakota")
SD_BHCF3 <- nless_mlatlongdf(df46, 44.506642, -103.914218, 44.445140, -103.624872, "South Dakota")
SD_BHCF4 <- greater_mlatlongdf(df46, 44.488946, -104.063920, 44.506642, -103.914218, "South Dakota")
SD_BHCF5 <- greater_mlatlongdf(df44_46, 43.536202, -103.503920, 44.129708, -103.303781, "South Dakota")
SD_BHCF6 <- rbind(SD_BHCF1, SD_BHCF2,SD_BHCF3, SD_BHCF5)
SD_BHCF7 <- nless_mlatlongdf(df44, 43.805954, -104.068539, 43.536202, -103.503920,  "South Dakota")
SD_BHCF8 <- anti_join(SD_BHCF6, SD_BHCF7)
SD_BHCF <- anti_join(SD_BHCF8, SD_BHCF4)

WY_BHCF1 <- ngreater_mlatlongdf(df46, 44.317627, -104.347260, 44.026411, -104.045762, "Wyoming")
WY_BHCF2 <- ngreater_mlatlongdf(df46, 44.507149, -104.118121, 44.026411, -104.045762, "Wyoming")
WY_BHCF3 <- ngreater_mlatlongdf(df46, 44.486114, -104.076571, 44.026411, -104.045762, "Wyoming")
WY_BHCF4 <- rbind(WY_BHCF1, WY_BHCF2)
WY_BHCF5 <- anti_join(WY_BHCF4, WY_BHCF3)
WY_BHCF6 <- less_mlatlongdf(df44, 44.486114, -104.076571, 44.488935, -104.052703, "Wyoming")
WY_BHCF <- rbind(WY_BHCF6, WY_BHCF5)

BHCF <- rbind(WY_BHCF, SD_BHCF)

## Southeastern Mixed Forest

AR_SMF <- LatLongdf(df36, 35.478800, -94.440780, 35.268162, -92.692572)
AR_SMF_shannon <- gen_shannon(AR_SMF)
OK_SMF <- LatLongdf(df36, 35.493475, -94.903364, 35.248540, -94.428764)
OK_SMF_shannon <- gen_shannon(OK_SMF)
SMF <- rbind(OK_SMF, AR_SMF)
SMF<-undup(SMF)

## Great Plains Steppe Province

ND_GPSP1 <- ngreater_mlatlongdf(df46_48, 49.100620, -103.209502, 45.836969, -98.201280, "North Dakota")
ND_GPSP <- long_less(ND_GPSP1, -97.816674, "North Dakota")

SD_GPSP1 <- less_mlatlongdf(df44_48, 42.993676, -100.169558, 46.119964, -98.382271, "South Dakota")
SD_GPSP <- long_less(SD_GPSP1, -97.771426, "South Dakota")

NE_GPSP1 <- less_mlatlongdf(df44, 42.294743, -101.730606, 43.109402, -99.649208, "Nebraska")
NE_GPSP2 <- ngreater_mlatlongdf(df40_42, 42.294743, -101.730606, 39.943141, -99.043452, "Nebraska")
NE_GPSP3 <- rbind(NE_GPSP1, NE_GPSP2)
NE_GPSP <- long_less(NE_GPSP3, -98.335079, "Nebraska")

KS_GPSP1 <- long_greater(df38_40, -99.455443, "Kansas")
KS_GPSP2 <- long_greater(KS_GPSP1, -97.304345, "Kansas")
KS_GPSP <- anti_join(KS_GPSP1, KS_GPSP2)

GPSP <- rbind(KS_GPSP, SD_GPSP, ND_GPSP, NE_GPSP)

## Northern Rocky Mountains Forest Steppe

WA_NRMFS1 <- ngreater_mlatlongdf(df48, 48.258087, -119.071937, 47.934973, -118.175979, "Washington")
WA_NRMFS2 <- ngreater_mlatlongdf(df48, 47.934973, -118.175979, 47.780591, -116.987194, "Washington")
WA_NRMFS3 <- ngreater_mlatlongdf(df50, 49.085789, -119.216958, 48.258087, -119.071937, "Washington")
WA_NRMFS <- rbind(WA_NRMFS1, WA_NRMFS2, WA_NRMFS3)

ID_NRMFS1<-ngreater_mlatlongdf(df48, 47.837628, -117.129509, 46.433551, -115.893530, "Idaho")
ID_NRMFS2 <- less_mlatlongdf(df48, 46.433551, -115.893530, 46.903083, -114.708337, "Idaho")
ID_NRMFS3 <- anti_join(ID_NRMFS1, ID_NRMFS2)
ID_NRMFS4 <- lat_greater(df50, 47.837628, "Idaho")
ID_NRMFS <- rbind(ID_NRMFS3, ID_NRMFS4)


MT_NRMFS1 <- greater_mlatlongdf(df48, 46.819537, -114.891614, 48.059089, -113.266556, "Montana")
MT_NRMFS2 <- nless_mlatlongdf(df50, 49.157757, -114.378438, 48.059089, -113.266556,  "Montana")
MT_NRMFS <- rbind (MT_NRMFS1, MT_NRMFS2)

NRMFS <- rbind(WA_NRMFS, ID_NRMFS, MT_NRMFS)


## Cascade Mixed Forest
df48_50 <- rbind(df48, df50)
WA_CMF1 <- nless_mlatlongdf(df48_50, 48.132279, -123.657060, 46.224089, -123.056276, "Washington")
WA_CMF2 <- ngreater_mlatlongdf(df50, 49.016499, -121.634011, 48.368716, -121.243846, "Washington")
WA_CMF3 <-less_mlatlongdf(df48_50, 47.551097, -121.291619, 48.368716, -121.243846, "Washington")
WA_CMF4 <- less_mlatlongdf(df46_48, 45.528304, -122.276301, 47.551097, -121.291619, "Washington")
WA_CMF5 <- rbind(WA_CMF1, WA_CMF2, WA_CMF3, WA_CMF4)
WA_CMF6 <- less_mlatlongdf(WA_CMF5, 45.595960, -121.450848, 48.058261, -120.132843, "Washington")
WA_CMF7 <- anti_join(WA_CMF5, WA_CMF6)
WA_CMF8 <- ngreater_mlatlongdf(WA_CMF5, 49.037733, -120.916186, 48.058261, -120.132843, "Washington")
WA_CMF <- anti_join(WA_CMF7, WA_CMF8)

OR_CMF1 <- greater_mlatlongdf(df46_48, 45.094673, -123.623736, 46.209311, -123.043191, "Oregon")
OR_CMF2 <- nless_mlatlongdf(df44_46, 45.094673, -123.623736, 43.889038, -123.212210, "Oregon")
OR_CMF3 <- less_mlatlongdf(df44_46, 43.889038, -123.212210, 45.343140, -122.021726, "Oregon")
OR_CMF4 <- ngreater_mlatlongdf(df46, 45.703541, -122.300975, 45.343140, -122.021726, "Oregon")
OR_CMF8 <- rbind(OR_CMF1, OR_CMF2, OR_CMF3, OR_CMF4)
OR_CMF9 <- ngreater_mlatlongdf(df46, 44.989036, -121.746691, 44.021293, -121.742476, "Oregon")
OR_CMF10 <- ngreater_mlatlongdf(df44_46, 45.770520, -121.846193, 44.989036, -121.746691, "Oregon")
OR_CMF11 <- anti_join(OR_CMF8, OR_CMF9)
OR_CMF11<- anti_join(OR_CMF11, OR_CMF10)
OR_CMF11 <- long_less(OR_CMF11, -121.742476, "Oregon")
OR_CMF12 <- lat_greater(OR_CMF11, 44.008377, "Oregon")
OR_CMF13 <- greater_mlatlongdf(df44_46, 43.552817, -124.257357,44.030218, -123.114766, "Oregon")
OR_CMF14 <- ngreater_mlatlongdf(df44_46, 44.030218, -123.114766, 43.622641, -122.815252, "Oregon")
OR_CMF15 <- less_mlatlongdf(df44, 43.115182, -123.094482, 43.622641, -122.815252, "Oregon")
OR_CMF16 <- ngreater_mlatlongdf(df44, 43.622641, -122.815252, 42.874784, -122.757080, "Oregon")
OR_CMF17 <- rbind(OR_CMF12, OR_CMF13, OR_CMF14, OR_CMF15, OR_CMF16)
OR_CMF18 <- less_mlatlongdf(OR_CMF17, 42.874784, -122.757080, 43.354640, -122.198620, "Oregon")
OR_CMF19 <- anti_join(OR_CMF17, OR_CMF18)
OR_CMF20 <- less_mlatlongdf(OR_CMF19, 43.354640, -122.198620, 44.006725, -121.846674, "Oregon")
OR_CMF21 <- anti_join(OR_CMF19, OR_CMF20)
OR_CMF <- long_less(OR_CMF21, -121.742476, "Oregon")

CAMF <- rbind(WA_CMF, OR_CMF)

## Pacific Lowland Mixed Forest
df46_50 <- rbind (df46, df48, df50)
WA_PLMF1 <- nless_mlatlongdf(df46_50, 48.485497, -123.028914, 45.959291, -122.851740, "Washington")

WA_PLMF2 <- nless_mlatlongdf(df46_48, 46.355175, -122.674567, 45.588563, -122.485963, "Washington") 
WA_PLMF3 <- greater_mlatlongdf(df46_50, 46.355175, -122.674567, 49.001818, -121.960159, "Washington")
WA_PLMF4 <- rbind(WA_PLMF2, WA_PLMF3)
WA_PLMF <- anti_join(WA_PLMF4, WA_PLMF1)


OR_PLMF1 <- nless_mlatlongdf(df46, 45.573204, -122.956450, 45.437981, -122.759864, "Oregon")
OR_PLMF2 <- less_mlatlongdf(OR_PLMF1, 45.412255, -123.089500, 45.948213, -122.799969, "Oregon")

OR_PLMF3 <- less_mlatlongdf(df46, 44.117834, -123.146809, 45.412255, -123.089500, "Oregon")
OR_PLMF4 <- greater_mlatlongdf(OR_PLMF3, 44.117834, -123.146809, 45.129567, -122.656791, "Oregon")

OR_PLMF5 <- ngreater_mlatlongdf(df46, 45.948213, -122.799969, 45.464268, -122.636577, "Oregon")
OR_PLMF6 <- nless_mlatlongdf(OR_PLMF5, 46.126857, -122.687076, 45.428483, -122.282759, "Oregon")

OR_PLMF8 <- less_mlatlongdf(df46, 45.161322, -123.196076, 45.527068, -123.128619, "Oregon")
OR_PLMF9 <- less_mlatlongdf(df46, 45.161322, -122.704605, 45.500556, -122.328774, "Oregon")
OR_PLMF10 <- anti_join(OR_PLMF8, OR_PLMF9)
OR_PLMF10 <- long_less(OR_PLMF10, -122.328774, "Oregon")


OR_PLMF11 <- rbind(OR_PLMF2, OR_PLMF4, OR_PLMF6, OR_PLMF10)
OR_PLMF <- undup(OR_PLMF11)

PLMF <- rbind(WA_PLMF, OR_PLMF)

## California COastal Chapparel Forest

OR_CCCF1 <- LatLongdf(df44, 43.448858, -124.299454, 42.020533, -122.568503)
OR_CCCF2 <- LatLongdf(df44, 42.622162, -124.438782, 42.009673, -124.046097)
OR_CCCF <- anti_join(OR_CCCF1, OR_CCCF2)


CA_CCCF1 <- less_mlatlongdf(df42_44, 40.659026, -123.814412, 42.042447, -123.304963, "California")
CA_CCCF2 <- ngreater_mlatlongdf(df40_42, 40.659026, -123.814412, 38.386713, -122.523894,  "California")
CA_CCCF3 <- rbind (CA_CCCF1, CA_CCCF2)
CA_CCCF4 <- nless_mlatlongdf(CA_CCCF3, 38.638817, -122.148003, 38.386713, -122.523894, "California")
CA_CCCF5 <- anti_join(CA_CCCF3, CA_CCCF4)
CA_CCCF6 <- ngreater_mlatlongdf(df40_42, 40.320744, -122.968703, 38.386713, -122.523894, "California")
CA_CCCF7 <- anti_join(CA_CCCF5, CA_CCCF6)
CA_CCCF8 <- less_mlatlongdf(df42, 40.320744, -122.968703, 40.538332, -122.536849,"California")
CA_CCCF9 <- anti_join(CA_CCCF7, CA_CCCF8)
CA_CCCF10 <- ngreater_mlatlongdf(df42, 40.538332, -122.536849, 40.513695, -121.815294, "California")
CA_CCCF11 <- ngreater_mlatlongdf(df40_42, 40.513695, -121.815294, 38.355219, -120.657615, "California")
CA_CCCF12 <- rbind(CA_CCCF10, CA_CCCF11, CA_CCCF9)
CA_CCCF13 <- ngreater_mlatlongdf(df38_40, 38.355219, -120.657615, 36.325875, -118.849188, "California")
CA_CCCF14 <- less_mlatlongdf(df36_38, 35.677253, -118.783928, 36.325875, -118.849188, "California")
CA_CCCF15 <- rbind(CA_CCCF12, CA_CCCF13, CA_CCCF14)

CA_CCCF16 <- less_mlatlongdf(df36_38, 35.677253, -118.783928, 36.152930, -118.140650, "California")
CA_CCCF17 <- ngreater_mlatlongdf(df38_40, 38.860033, -120.110588, 36.152930, -118.140650, "California")
CA_CCCF20 <- rbind(CA_CCCF16, CA_CCCF17)
CA_CCCF21 <- anti_join(CA_CCCF15, CA_CCCF20)

CA_CCCF22 <- ngreater_mlatlongdf(df40_42, 40.427168, -120.901869, 39.889832, -119.967635, "California")
CA_CCCF23 <- less_mlatlongdf(df42, 40.427168, -120.901869, 41.011138, -119.967635, "California")
CA_CCCF24 <- rbind(CA_CCCF22, CA_CCCF23)
CA_CCCF <- anti_join(CA_CCCF21, CA_CCCF24)

CCCF <- rbind(OR_CCCF, CA_CCCF)


## California Coastal Range Open Woodland

CA_CCROW1 <- ngreater_mlatlongdf(df36_38, 37.433612, -121.725313, 35.194035, -120.087835, "California")
CA_CCROW2 <- ngreater_mlatlongdf(df36_38, 37.433612, -121.794652, 34.857433, -118.965451, "California")
CA_CCROW3 <- anti_join(CA_CCROW1, CA_CCROW2)

CA_CCROW4 <- ngreater_mlatlongdf(df36, 35.194035, -120.087835, 34.460152, -118.819905, "California")
CA_CCROW5 <- ngreater_mlatlongdf(df36, 34.460152, -118.819905, 34.208488, -118.042966, "California")
CA_CCROW6 <- ngreater_mlatlongdf(df34_36, 34.208488, -118.042966, 33.254225, -117.095066, "California")
CA_CCROW7 <- ngreater_mlatlongdf(df34, 33.254225, -117.095066, 32.546944, -116.568634, "California")
CA_CCROW8 <- rbind(CA_CCROW4, CA_CCROW5, CA_CCROW6, CA_CCROW7)
CA_CCROW13 <- anti_join(CA_CCROW3, CA_CCROW8)

CA_CCROW9 <- ngreater_mlatlongdf(df34, 33.127023, -116.630593, 32.651339, -116.081815, "California")
CA_CCROW10 <- ngreater_mlatlongdf(df34_36,34.217195, -116.745660, 33.127023, -116.630593, "California")
CA_CCROW11 <- ngreater_mlatlongdf(df36, 34.938671, -119.321376, 34.217195, -116.745660, "California")
CA_CCROW14 <- ngreater_mlatlongdf(df36, 35.194035, -119.746237, 34.938671, -119.321376, "California")
CA_CCROW12 <- rbind (CA_CCROW11, CA_CCROW10, CA_CCROW9, CA_CCROW14)

CA_CCROW13 <- anti_join(CA_CCROW8, CA_CCROW12)
CA_CCROW <- rbind(CA_CCROW13, CA_CCROW3)

## Middle Rocky Mountain Steppe

ID_MRMS1 <- LatLongdf(df46, 45.272317, -116.093996, 44.228952, -113.225519)
ID_MRMS <- ID_MRMS1[ID_MRMS1$stateProvince == "Idaho", ]

MT_MRMS1 <- LatLongdf(df46_48, 46.538008, -114.309943, 45.260006, -111.826261)
MT_MRMS <- MT_MRMS1[MT_MRMS1$stateProvince == "Montana", ]

OR_MRMS1 <- less_mlatlongdf(df46, 44.553911, -120.868961, 45.444379, -118.490223, "Oregon")
OR_MRMS2 <- less_mlatlongdf(df46_48, 45.444379, -118.490223, 46.212188, -117.825576, "Oregon")
OR_MRMS3 <- less_mlatlongdf(df46, 44.241484, -119.259815, 45.444379, -118.490223,  "Oregon")
OR_MRMS4 <- rbind(OR_MRMS1, OR_MRMS2)

OR_MRMS5 <- less_mlatlongdf(df46, 44.241484, -119.259815, 45.333828, -116.461301, "Oregon")
OR_MRMS6 <- anti_join(OR_MRMS3, OR_MRMS5)
OR_MRMS <- rbind(OR_MRMS4, OR_MRMS6)

WA_MRMS <- less_mlatlongdf(df46_48, 45.444379, -118.490223, 46.212188, -117.825576, "Washington")

MRMS <- rbind(ID_MRMS, MT_MRMS, OR_MRMS, WA_MRMS)


## Great Plains Palousse Dry Steppe

WA_GPPS <- LatLongdf(df48, 47.632617, -119.749026, 46.424603, -117.863897)

OR_GPPS1 <- greater_mlatlongdf(df46, 45.018630, -120.538524, 45.517568, -118.975640, "Oregon")
OR_GPPS2 <- greater_mlatlongdf(df46_48, 45.517568, -118.975640, 46.056864, -118.234478, "Oregon")
OR_GPPS3 <- greater_mlatlongdf(df46, 44.584199, -121.070227, 45.018630, -120.538524, "Oregon")
OR_GPPS4 <- rbind(OR_GPPS1, OR_GPPS2, OR_GPPS3)

OR_GPPS5 <- LatLongdf(df44, 43.577266, -120.490187, 42.012734, -117.037422)

OR_GPPS6 <- less_mlatlongdf(df44_46, 43.577266, -118.548870, 44.451885, -117.195005, "Oregon")
OR_GPPS7 <- rbind(OR_GPPS6, OR_GPPS5, OR_GPPS4)                            
OR_GPPS <- long_greater(OR_GPPS7, -121.171467, "Oregon")
                    


NV_GPPS1 <- greater_mlatlongdf(df42_44, 40.424779, -120.086540, 42.114975, -118.563799, "Nevada")
NV_GPPS2 <- LatLongdf(df42_44, 42.097594, -117.345606, 41.573969, -113.972149)
NV_GPPS3 <- NV_GPPS2[NV_GPPS2$stateProvince == "Nevada", ]

NV_GPPS <- rbind(NV_GPPS3, NV_GPPS1)

ID_GPPS1 <- nless_mlatlongdf(df44_46, 44.392384, -117.251899, 43.200322, -115.635451, "Idaho")
ID_GPPS2 <- less_mlatlongdf(df44, 43.260063, -113.269346, 43.955527, -112.121433, "Idaho")
ID_GPPS3 <- lat_less(df42_44, 43.200322, "Idaho")
ID_GPPS4 <- rbind(ID_GPPS1, ID_GPPS2, ID_GPPS3)

ID_GPPS5 <- less_mlatlongdf(df42_44, 41.993212, -114.663239, 43.710489, -111.746605, "Idaho")
ID_GPPS6 <- anti_join(ID_GPPS4, ID_GPPS5)

ID_GPPS <- long_less(ID_GPPS6, -111.559190, "Idaho")


WY_GPPS1 <- nless_mlatlongdf(df44_46, 45.168336, -108.409854, 43.987674, -107.322034, "Wyoming")
WY_GPPS2 <- greater_mlatlongdf(df44, 43.369531, -107.969545, 43.987674, -107.322034, "Wyoming")
WY_GPPS3 <- nless_mlatlongdf(df44, 43.369531, -107.969545, 42.459054, -106.311915, "Wyoming")

WY_GPPS4 <- rbind(WY_GPPS1, WY_GPPS2, WY_GPPS3)

wy_test <- draw_complex(43.987674, -107.322034, 45.168336, -108.409854, "nless", 43.987674, -107.322034, 43.369531, -107.969545, "greater", 42.459054, -106.311915, 43.369531, -107.969545, "nless", "Wyoming")

WY_GPPS5 <- draw_complex(42.459054, -106.311915, 41.633456, -107.510138, "greater", 41.633456, -107.510138, 40.990010, -107.990382, "greater", 40.990010, -107.990382, 40.960611, -108.885971, "greater", "Wyoming")

WY_GPPS6 <- rbind(WY_GPPS5, wy_test)

WY_GPPS7 <- draw_complex(45.168336, -108.409854, 43.478507, -108.989808, "greater", 42.424029, -108.600421, 43.478507, -108.989808, "nless", 40.960611, -108.885971, 42.682186, -109.911356, "nless", "Wyoming")

WY_GPPS8 <- anti_join(WY_GPPS6, WY_GPPS7)
GPPS <- rbind(WY_GPPS8, ID_GPPS, NV_GPPS, OR_GPPS, WA_GPPS)

## Southern Rocky Mountain Steppe-Open Woodland

WY_SRMSOW1 <- LatLongdf(df46, 45.413600, -111.251831,44.178067, -109.745917)
WY_SRMSOW2 <- greater_mlatlongdf(df42_44, 42.558932, -111.204771, 44.178067, -109.745917, "Wyoming")
WY_SRMSOW3 <- less_mlatlongdf(df42, 40.965392, -107.234650, 41.611880, -106.538237, "Wyoming")
WY_SRMSOW4 <- nless_mlatlongdf(WY_SRMSOW3, 41.611880, -106.538237, 40.984240, -106.133440, "Wyoming")
WY_SRMSOW <- rbind(WY_SRMSOW1, WY_SRMSOW2, WY_SRMSOW4)

CO_SRMSOW1 <- draw_complex(39.767415, -105.356127, 41.021929, -105.495852, "nless", 38.719496, -104.900180, 39.767415, -105.356127, "nlessr", 38.719496, -104.900180, 37.885808, -105.027168, "greater", "Colorado")

CO_SRMSOW2 <- nless_mlatlongdf(df38, 37.885808, -105.027168, 37.000496, -104.976939, "Colorado")
CO_SRMSOW3 <- rbind(CO_SRMSOW1, CO_SRMSOW2)

CO_SRMSOW4 <- draw_complex(41.041424, -107.518534, 39.763942, -108.071054, "greater", 39.763942, -108.071054, 38.963981, -108.342292, "greater", 38.254634, -107.696239, 38.963981, -108.342292, "nless", "Colorado")
CO_SRMSOW5 <- draw_complex(38.254634, -107.696239, 37.853122, -108.196498, "greater", 37.300602, -107.376309, 37.853122, -108.196498, "nless", 36.984598, -107.205436, 37.300602, -107.376309, "nless", "Colorado")
CO_SRMSOW7 <- rbind(CO_SRMSOW5, CO_SRMSOW4)

CO_SRMSOW12 <- anti_join(CO_SRMSOW3, CO_SRMSOW7)

CO_SRMSOW6 <- ngreater_mlatlongdf(df38, 37.720396, -106.299428, 36.995594, -106.071459, "Colorado")
CO_SRMSOW8 <- less_mlatlongdf(df38_40, 37.720396, -106.299428, 38.348838, -106.109454, "Colorado")
CO_SRMSOW9 <- rbind(CO_SRMSOW6, CO_SRMSOW8)

CO_SRMSOW10 <- nless_mlatlongdf(CO_SRMSOW9, 38.348838, -106.109454, 36.965242, -105.438213, "Colorado")

CO_SRMSOW <- anti_join(CO_SRMSOW12, CO_SRMSOW10)                             

NM_SRMSOW1 <- less_mlatlongdf(df36_38, 35.926103, -105.970140, 37.056261, -105.488873, "New Mexico")
NM_SRMSOW2 <- greater_mlatlongdf(NM_SRMSOW1, 35.926103, -105.970140, 36.304654, -105.362224, "New Mexico")
NM_SRMSOW3 <- greater_mlatlongdf(NM_SRMSOW1, 36.304654, -105.362224, 37.005708, -105.349559, "New Mexico")
NM_SRMSOW4 <- rbind(NM_SRMSOW2, NM_SRMSOW3)

NM_SRMSOW5 <- nless_mlatlongdf(df38, 37.066368, -106.350088, 36.528871, -106.033465, "New Mexico")
NM_SRMSOW6 <- greater_mlatlongdf(df34_36, 35.751564, -106.818690, 36.528871, -106.033465, "New Mexico")
NM_SRMSOW7 <- rbind(NM_SRMSOW5, NM_SRMSOW6)
NM_SRMSOW8 <- ngreater_mlatlongdf(NM_SRMSOW7, 37.076473, -106.983334, 35.751564, -106.818690, "New Mexico")
NM_SRMSOW <- rbind(NM_SRMSOW8, NM_SRMSOW4)

UT_SRMSOW1 <- draw_complex(40.186533, -111.599060, 42.079519, -112.174311, "ngreater", 40.186533, -111.599060, 39.918507, -111.761631, "less", 39.707176, -111.136357, 39.918507, -111.761631, "ngreater", "Utah")

UT_SRMSOW3 <- draw_complex(40.955948, -111.186379, 42.051668, -111.549038, "ngreater", 40.851979, -109.898316, 40.955948, -111.186379, "ngreater", 40.851979, -109.898316, 40.415430, -110.986292, "less", "Utah")
UT_SRMSOW4 <- less_mlatlongdf(df40_42, 39.707176, -111.136357, 40.415430, -110.986292, "Utah")
UT_SRMSOW5 <- rbind(UT_SRMSOW3, UT_SRMSOW4)
UT_SRMSOW <- anti_join(UT_SRMSOW1, UT_SRMSOW5)

ID_SRMSOW1 <- less_mlatlongdf(df42_44, 41.968040, -112.699541, 43.510928, -111.173874, "Idaho")


SRMSOW <- rbind(UT_SRMSOW, NM_SRMSOW, CO_SRMSOW, WY_SRMSOW, ID_SRMSOW1)

## California Valley Grassland

CVG1 <- draw_complex(39.939424, -122.178413,39.040442, -122.066402,"less",36.957600, -120.789480, 39.040442, -122.066402, "ngreater", 35.375233, -119.400548, 36.957600, -120.789480, "ngreater", "California")
CVG2 <- draw_complex(39.179503, -121.506349, 39.939424, -122.178413, "ngreater", 39.179503, -121.506349, 38.481458, -121.584756, "less", 36.346535, -119.355744, 38.481458, -121.584756, "ngreater", "California")
CVG3 <- anti_join(CVG1, CVG2)                     
CVG4 <-less_mlatlongdf(df36_38, 35.375233, -119.400548, 36.346535, -119.355744, "California")
CVG<- anti_join(CVG3, CVG4)

## American Semi-Desert
CA_ASD1 <- draw_complex(36.818736, -116.910322, 36.306429, -117.833750, "less", 36.306429, -117.833750, 34.725789, -117.843270, "less", 33.578147, -114.474727, 34.725789, -117.843270, "ngreater", "California")
CA_ASD2 <- ngreater_mlatlongdf(df34_36, 34.418847, -116.409461, 32.658158, -115.483761, "California")
CA_ASD3 <- rbind(CA_ASD1, CA_ASD2)
CA_ASD <- undup(CA_ASD3)

NV_ASD1 <- draw_complex(36.699923, -116.381177, 36.929384, -117.239888, "nless", 36.358068, -115.346347, 36.929384, -117.239888, "nless", 36.098363, -114.713951, 36.358068, -115.346347, "nless", "Nevada")
NV_ASD2 <- lat_less(df36_38, 36.098363, "Nevada")
NV_ASD3 <- less_mlatlongdf(df38, 36.098363, -114.713951, 37.296812, -114.024065, "Nevada")
NV_ASD <- rbind(NV_ASD1, NV_ASD2, NV_ASD3)

AZ_ASD <- draw_complex(34.251217, -113.499355, 36.198000, -114.442683, "nless", 33.196579, -111.510257, 34.251217, -113.499355, "nless", 33.196579, -111.510257, 31.445748, -111.568273,  "greater", "Arizona")

ASD <- rbind(CA_ASD, NV_ASD, AZ_ASD)

##Arizona-New Mexico Semi-Desert Mountains

AZ_ADM1 <- draw_complex(34.078550, -111.2310, 35.518859, -111.990119, "ngreater", 33.267276, -109.750283, 34.078550, -111.2310, "ngreater", 33.215673, -108.997550, 33.267276, -109.750283, "ngreater", "Arizona")
AZ_ADM2 <- draw_complex(34.664199, -111.218729, 35.518859, -111.990119, "ngreater", 34.196007, -110.213029, 34.664199, -111.218729, "ngreater", 33.996745, -109.028399, 34.196007, -110.213029, "ngreater", "Arizona")
AZ_ADM <- anti_join(AZ_ADM1, AZ_ADM2)

NM_ADM1 <- draw_complex(35.660263, -108.521701, 36.449273, -109.090373, "nless", 35.228282, -107.445835, 35.660263, -108.521701, "nless",33.975893, -107.084652, 35.228282, -107.445835, "nless", "New Mexico")
NM_ADM2 <- draw_complex(33.975893, -107.084652, 32.963064, -107.684063, "greater", 32.963064, -107.684063, 32.737108, -108.490962, "greater", 32.737108, -108.490962, 32.620678, -109.105742, "greater", "new Mexico")
NM_ADM3 <- rbind(NM_ADM1, NM_ADM2)

NM_ADM4 <- draw_complex(34.806610, -108.544755, 35.203169, -109.128797, "nless", 34.806610, -108.544755, 34.560160, -108.575494, "greater", 34.560160, -108.575494, 34.109616, -109.090373, "greater", "New Mexico")
NM_ADM5 <- anti_join(NM_ADM3, NM_ADM4)

NM_ADM6 <- draw_complex(36.190202, -104.469010, 34.883343, -105.533576, "less", 34.883343, -105.533576,33.408461, -106.688234, "less", 31.977333, -106.405845, 33.408461, -106.688234, "ngreater", "New Mexico")
NM_ADM7 <- draw_complex(36.190202, -104.469010, 35.729774, -104.11465, "ngreater", 35.729774, -104.11465, 34.946263, -105.737230, "less", 32.825383, -105.532077, 34.946263, -105.737230, "ngreater", "New Mexico")
NM_ADM8 <- anti_join(NM_ADM6, NM_ADM7)
NM_ADM9 <- ngreater_mlatlongdf(NM_ADM8, 34.946263, -105.737230, 31.998790, -104.758089, "New Mexico")
NM_ADM10 <- anti_join(NM_ADM8, NM_ADM9)

NM_ADM11<-rbind(NM_ADM10, NM_ADM5)
NM_ADM <- long_less(NM_ADM11, -104.114654, "New Mexico")

ADM <- rbind(AZ_ADM, NM_ADM)

## Chihuahuan Semi_Desert

AZ_CSD <- draw_complex(32.702546, -111.265197, 31.425002, -111.498510, "less", 33.229555, -110.314000, 32.702546, -111.265197, "less", 33.229555, -110.314000, 32.551396, -108.967966, "nless", "Arizona")

NM_CSD1 <- nless_mlatlongdf(df34, 32.050794, -107.927033, 32.460583, -109.111543,"New Mexico")
NM_CSD2 <- nless_mlatlongdf(df32_34, 32.050794, -107.927033, 31.989928, -107.460408, "New Mexico")
NM_CSD3 <- rbind(NM_CSD1, NM_CSD2)
NM_CSD4 <- lat_less(df32, 31.989928, "New Mexico")
NM_CSD5 <- rbind(NM_CSD3, NM_CSD4)

NM_CSD6 <- less_mlatlongdf(df32_34, 31.989928, -107.460408,32.364224, -107.011396, "New Mexico")
NM_CSD7 <- ngreater_mlatlongdf(df34, 32.364224, -107.011396, 32.125650, -106.635293,"New Mexico")
NM_CSD8 <- less_mlatlongdf(df32_34, 31.779940, -106.748124, 32.125650, -106.635293, "New Mexico")
NM_CSD9 <- rbind(NM_CSD8, NM_CSD7)
NM_CSD10 <- anti_join(NM_CSD9, NM_CSD6)

NM_CSD11 <- draw_complex(32.201216, -104.41477, 31.947421, -104.437145, "less", 32.201216, -104.41477, 32.729670, -105.041121, "ngreater", 33.642300, -105.141784, 32.729670, -105.041121, "ngreater", "New Mexico")
NM_CSD12 <- draw_complex(33.259686, -104.588139, 33.642300, -105.141784, "ngreater", 33.259686, -104.588139, 32.739078, -104.616101, "less", 31.992763, -104.006532, 33.259686, -104.588139, "ngreater", "New Mexico")
NM_CSD13 <- anti_join(NM_CSD11, NM_CSD12)
 NM_CSD <- rbind(NM_CSD13, NM_CSD10, NM_CSD5)

TX_CSD1 <- long_less(df32_34, -103.099927, "Texas")
TX_CSD2 <- ngreater_mlatlongdf(df32_34, 32.049403, -105.298374, 30.993621, -104.571612, "Texas")
TX_CSD3 <- greater_mlatlongdf(TX_CSD2, 30.993621, -104.571612, 32.049403, -104.469251, "Texas")
TX_CSD4 <- anti_join(TX_CSD1, TX_CSD3)
TX_CSD5 <- nless_mlatlongdf(df32_34, 32.069674, -103.099927, 29.399220, -101.062824, "Texas")
TX_CSD <- rbind(TX_CSD5, TX_CSD4)

CSD <- rbind(TX_CSD, NM_CSD, AZ_CSD)

## Colorado Plateau Semi Arid Desert

AZ_CPAD1 <- draw_complex(32.993808, -111.352296, 34.408899, -113.089182, "ngreater", 34.408899, -113.089182, 35.561736, -113.6980, "ngreater", 35.561736, -113.6980, 36.372083, -114.076892, "ngreater", "Arizona")
AZ_CPAD2 <- draw_complex(32.993808, -111.352296, 34.408899, -112.182653, "ngreater", 34.408899, -112.182653, 35.197692, -112.534440, "ngreater", 35.583747, -111.939108, 35.197692, -112.534440, "less", "Arizona")
AZ_CPAD3 <- anti_join(AZ_CPAD1, AZ_CPAD2)

AZ_CPAD4 <- draw_complex(34.843114, -110.883747, 35.638746, -111.857927, "ngreater", 34.843114, -110.883747, 34.565033, -109.760734, "ngreater", 34.140563, -108.975978, 34.565033, -109.760734, "ngreater", "Arizona")

AZ_CPAD5 <- rbind(AZ_CPAD4, AZ_CPAD3)
AZ_CPAD6 <- lat_greater(df38, 36.372083, "Arizona")
AZ_CPAD7 <- rbind(AZ_CPAD5, AZ_CPAD6)
AZ_CPAD8 <- draw_complex(37.074197, -109.337790, 35.872722, -109.575101, "less", 35.872722, -109.575101, 35.551577, -109.733309, "less", 35.128565, -109.012585, 35.551577, -109.733309, "ngreater", "Arizona")
AZ_CPAD <- anti_join(AZ_CPAD7, AZ_CPAD8)

NM_CPAD1 <- draw_complex(35.902056, -108.564314, 37.104618, -109.241609, "ngreater", 35.493308, -107.637805, 35.902056, -108.564314, "ngreater", 35.493308, -107.637805, 34.931734, -107.062675, "ngreater", "New Mexico")
NM_CPAD2 <- draw_complex(35.931812, -107.315171, 37.037061, -107.553640, "ngreater", 35.590328, -106.971496, 35.931812, -107.315171, "ngreater", 35.590328, -106.971496, 34.931734, -107.062675, "less", "New Mexico")
NM_CPAD3 <- anti_join(NM_CPAD1, NM_CPAD2)
NM_CPAD4 <- less_mlatlongdf(df36, 35.441897, -107.665860, 35.635944, -107.22399, "New Mexico")
NM_CPAD5 <- ngreater_mlatlongdf(df36, 35.635944, -107.22399, 35.190082, -107.167882, "New Mexico")
NM_CPAD6 <- anti_join(NM_CPAD4, NM_CPAD5)
NM_CPAD7 <- anti_join(NM_CPAD3, NM_CPAD6)

CO_CPAD <- nless_mlatlongdf(df38, 37.737345, -109.079165, 36.982379, -108.282312, "Colorado")

UT_CPAD <- less_mlatlongdf(df38, 37.009474, -112.030662, 37.539193, -109.015456, "Utah")

CPAD <- rbind(UT_CPAD, CO_CPAD, NM_CPAD7, AZ_CPAD)

##Southwest Plateau

NM_SP <- LatLongdf(df34_36, 35.165392, -104.398591, 32.021491, -103.096342)
df32_38 <- rbind(df32_36, df38)
TX_SP1 <- ngreater_mlatlongdf(df32_34, 32.162225, -103.108647, 28.752924, -100.543045, "Texas")
TX_SP2 <- ngreater_mlatlongdf(TX_SP1, 34.114783, -100.644394, 29.389513, -97.678064, "Texas")
TX_SP3 <- LatLongdf(df32_38, 36.087745, -103.017932, 31.966697, -99.613679)
TX_SP5 <- anti_join(TX_SP1, TX_SP2)
TX_SP4 <- rbind(TX_SP5, TX_SP3)
TX_SP6 <- undup(TX_SP4)
TX_SP7 <- TX_SP6[TX_SP6$stateProvince == "Texas", ]
TX_SP8 <- less_mlatlongdf(df36, 34.817670, -100.73630, 35.217351, -99.971725, "Texas")
TX_SP9 <- ngreater_mlatlongdf(df34_36, 34.817670, -100.73630, 34.456291, -99.662641, "Texas")
TX_SP10 <- rbind(TX_SP8, TX_SP9)
TX_SP <- anti_join(TX_SP7, TX_SP10)

SP <- rbind(NM_SP, TX_SP)

## Subtropical Prairie Parkland 

TX_SPP1 <- draw_complex(34.063246, -97.644653, 32.367515, -98.850385, "less", 31.507435, -96.631139, 32.367515, -98.850385, "ngreater", 31.507435, -96.631139, 30.353285, -97.277691, "less", "Texas")
TX_SPP <- long_less(TX_SPP1, -95.938942, "Texas")

OK_SPP <- LatLongdf(df34_36, 35.360484, -97.124665, 33.908735, -96.070862)

SPP <- rbind(OK_SPP, TX_SPP)

## Nevada-Utah Mountains

NV_NVUTM1 <- draw_complex(40.010267, -115.661102,39.648764, -116.678150, "less", 39.648764, -116.678150, 38.878897, -117.517149, "less", 38.344720, -116.429558, 38.878897, -117.517149, "ngreater", "Nevada")
NV_NVUTM2 <- draw_complex(38.915174, -115.419652, 38.878897, -117.517149, "less", 38.915174, -115.419652, 39.011820, -114.642801, "less", 39.011820, -114.642801, 40.010267, -115.661102, "ngreater", "Nevada")
NV_NVUTM <- anti_join(NV_NVUTM1, NV_NVUTM2)

UT_NVUTM1 <- draw_complex(39.082189, -112.219110, 38.296019, -112.964175, "less", 38.296019, -112.964175, 37.916031, -114.067976, "less", 37.435534, -113.847215, 37.916031, -114.067976, "ngreater", "Utah")
UT_NVUTM2 <- draw_complex(37.534070, -112.702023, 37.435534, -113.847215, "less", 38.046536, -111.943160, 37.534070, -112.702023, "less", 38.046536, -111.943160,39.082189, -112.219110, "ngreater", "Utah")
UT_NVUTM3 <- anti_join(UT_NVUTM1, UT_NVUTM2)
UT_NVUTM4 <- less_mlatlongdf(df40, 39.092898, -109.997711, 39.732482, -109.004291, "Utah")
UT_NVUTM5 <- less_mlatlongdf(df40, 39.092898, -109.997711,39.306744, -108.990494, "Utah")
UT_NVUTM6 <- anti_join(UT_NVUTM4, UT_NVUTM5)
UT_NVUTM <- rbind(UT_NVUTM6, UT_NVUTM3)

CO_NVUTM <- LatLongdf(df40, 39.933792, -108.838721,39.413423, -108.259226)
NVUTM <- rbind(CO_NVUTM, UT_NVUTM, NV_NVUTM)

## Intermountain Semi Desert

NV_ISD1 <- draw_complex(41.172771, -115.445154,39.143324, -118.709203, "greater", 37.666848, -116.513388, 39.143324, -118.709203, "nless", 36.816455, -114.035679, 37.666848, -116.513388, "nless", "Nevada")
NV_ISD2 <- lat_less(df36_38, 36.816455, "Nevada")
NV_ISD3 <- rbind(NV_ISD2, NV_ISD1)
NV_ISD4 <- greater_mlatlongdf(df42_44, 40.023760, -120.059332, 42.027049, -118.041556,"Nevada")
NV_ISD <- anti_join(NV_ISD3, NV_ISD4)

UT_ISD1 <- greater_mlatlongdf(df40_42, 38.619141, -114.116902,40.750311, -112.220545,"Utah")
UT_ISD2 <- LatLongdf(df42, 41.548817, -114.042046, 40.750311, -112.220545)
UT_ISD3 <- LatLongdf(df40_42, 40.361999, -114.043124,39.672874, -113.738769)
UT_ISD4 <- rbind(UT_ISD1, UT_ISD2)
UT_ISD5 <- anti_join(UT_ISD4, UT_ISD3)

UT_ISD6 <- LatLongdf(df40, 38.847164, -110.914702,38.268872, -109.470033)
UT_ISD7 <- less_mlatlongdf(df40, 38.847164, -110.914702,39.473916, -110.691122, "Utah")
UT_ISD8 <- ngreater_mlatlongdf(df40, 39.473916, -110.691122, 38.847164, -110.278360, "Utah")
UT_ISD9 <- anti_join(UT_ISD7, UT_ISD8)
UT_ISD10 <- rbind(UT_ISD9, UT_ISD6)
UT_ISD11 <- draw_complex(40.442614, -109.645398, 40.200032, -110.212946, "less", 40.121169, -109.224036, 40.200032, -110.212946, "ngreater", 40.121169, -109.224036, 40.442614, -109.645398, "ngreater", "Utah")
UT_ISD11 <- undup(UT_ISD11)

UT_ISD <- rbind(UT_ISD5, UT_ISD10, UT_ISD11)

CO_ISD <- LatLongdf(df38, 37.335069, -109.043982, 37.007975, -108.734428)
ISD <- rbind(UT_ISD, NV_ISD, CO_ISD)

## Data Transformation

MRF_shannon <- gen_shannon(Mississippi_River_Forest)
MF_shannon <- gen_shannon(Mixed_Forest)
CMF_shannon <- gen_shannon(Coastal_Mixed_Forest)
ABF_shannon <- gen_shannon(Appalacian_Broadleaf_Forest)
A_shannon <- gen_shannon(Adirondack)
LMF_shannon <- gen_shannon(LMF)
EBF_shannon <- gen_shannon(EBF)
EBFC_shannon <- gen_shannon(EBFC)
TPP_shannon <- gen_shannon(TPP)
GPS_shannon <- gen_shannon(GPS)
BHCF_shannon <- gen_shannon(BHCF)
SMF_shannon <- gen_shannon(SMF)
GPSP_shannon <- gen_shannon(GPSP)
NRMFs_shannon <- gen_shannon(NRMFS)
CAMF_shannon <- gen_shannon(CAMF)
PLMF_shannon <- gen_shannon(PLMF)
CCCF_shannon <- gen_shannon(CCCF)
CA_CCROW_shannon <- gen_shannon(CA_CCROW)
MRMS_shannon <- gen_shannon(MRMS)
GPPS_shannon <- gen_shannon(GPPS)
SRMSOW_shannon <- gen_shannon(SRMSOW)
CVG_shannon <- gen_shannon(CVG)
ASD_shannon <- gen_shannon(ASD)
ADM_shannon <- gen_shannon(ADM)
CSD_shannon <- gen_shannon(CSD)
CPAD_shannon <- gen_shannon(CPAD)
SP_shannon <- gen_shannon(SP)
SPP_shannon <- gen_shannon(SPP)
NVUTM_shannon <- gen_shannon(NVUTM)
ISD_shannon <- gen_shannon(ISD)

## Analysis
gen_m0 <- function(a){
  m0 = lm(Shannon ~ biome + lat + biome*lat, data = a)
  return(summary(m0))
}
gen_m1 <- function(a){
  m0 = lm(Shannon ~ lat, data = a)
  return(summary(m0))
}
FM <- function(x, y){
  r1=filter(x, Shannon_Diversity_Index != "N/A")
  r1$Shannon_Diversity_Index <- as.numeric(r1$Shannon_Diversity_Index)
  if(grepl("latband", colnames(x)[2])){
    r2<- r1 %>%
      group_by(latband) %>%
      summarize(Shannon = mean(exp(Shannon_Diversity_Index)))
    r2$latband <- str_remove(r2$latband, "band")
    col_names <- colnames(r2)
    col_names[1] <- "lat"
    colnames(r2) <- col_names
  }else{
    r2<- r1 %>%
      group_by(lat) %>%
      summarize(Shannon = mean(exp(Shannon_Diversity_Index)))
  }
  r3 <- r2 %>%
    mutate(r2, biome = rep(c(y)))
  return(as.data.frame(r3))
}

###Eastern Forests
TN_EBFC_shannon <- gen_shannon(TN_EBFC)
KY_EBFC_shannon <- gen_shannon(KY_EBFC)
OH_EBFC_shannon <- gen_shannon(OH_EBFC)
IN_EBFC1_shannon <- gen_shannon(IN_EBFC1)
IN_EBFC2_shannon <- gen_shannon(IN_EBFC2)
IL_EBFC1_shannon <- gen_shannon(IL_EBFC1)
IL_EBFC2_shannon <- gen_shannon(IL_EBFC2)
MO_EBFC_shannon <- gen_shannon(MO_EBFC)
IA_EBFC_shannon <- gen_shannon(IA_EBFC)
MN_EBFC_shannon <- gen_shannon(MN_EBFC)
WI_EBFC_shannon <- gen_shannon(WI_EBFC)
MI_EBFC_shannon <- gen_shannon(MI_EBFC)
AR_EBFC_shannon <- gen_shannon(AR_EBFC)

EBFC_shannon2 <- bind_rows(FM(TN_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(KY_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(OH_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(IN_EBFC1_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(IN_EBFC2_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(IL_EBFC1_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(IL_EBFC2_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(MO_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(IA_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(MN_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(WI_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(MI_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                           FM(AR_EBFC_shannon, "Eastern Broadleaf Forest (Continental)"))

EBF_shannon2 <- bind_rows(FM(PA_EBF_shannon, "Eastern Broadleaf Forest"), FM(MA_EBF_shannon, "Eastern Broadleaf Forest"), FM(TN_EBF_shannon, "Eastern Broadleaf Forest"), FM(KY_EBF_shannon, "Eastern Broadleaf Forest"), FM(WV_EBF_shannon, "Eastern Broadleaf Forest"), FM(OH_EBF_shannon, "Eastern Broadleaf Forest"), FM(NH_EBF_shannon, "Eastern Broadleaf Forest"), FM(ME_EBF_shannon, "Eastern Broadleaf Forest"))
MRF_shannon2 <- bind_rows(FM(LA_MRF_shannon, "Mississippi River Forest"), FM(MI_MRF_shannon, "Mississippi River Forest"), FM(AK_MRF_shannon, "Mississippi River Forest"), FM(MO_MRF_shannon, "Mississippi River Forest"), FM(TN_MRF_shannon,"Mississippi River Forest"), FM(KY_MRF_Shannon, "Mississippi River Forest"))
MF_shannon2 <- bind_rows(FM(AL_MF_shannon, "Mixed Forest"), FM(AK_MF_shannon,  "Mixed Forest"), FM(LA_MF_shannon, "Mixed Forest"), FM(TX_MF_shannon,  "Mixed Forest"), FM(GA_MF_shannon,  "Mixed Forest"), FM(NC_MF_shannon, "Mixed Forest"), FM(VA_MF_shannon,  "Mixed Forest"), FM(SC_MF_shannon,  "Mixed Forest"), FM(MI_MF_shannon,  "Mixed Forest"))
CMF_shannon2 <- bind_rows(FM(DA_CMF_shannon, "Coastal Mixed Forest"), FM(MD_CMF_shannon, "Coastal Mixed Forest"), FM(GA_CMF_shannon, "Coastal Mixed Forest"), FM(VA_CMF_shannon, "Coastal Mixed Forest"), FM(SC_CMF_shannon, "Coastal Mixed Forest"), FM(NC_CMF_shannon, "Coastal Mixed Forest"), FM(AL_CMF_shannon, "Coastal Mixed Forest"), FM(MI_CMF_shannon, "Coastal Mixed Forest"), FM(LA_CMF_shannon, "Coastal Mixed Forest"))
SMF_shannon2 <- bind_rows(FM(OK_SMF_shannon, "Southeast Mixed Forest"), FM(AR_SMF_shannon, "Southeast Mixed Forest"))
LMF_shannon2 <- bind_rows(FM(PA_LMF_shannon, "Laurentian_Mixed_Forest"), FM(NY_LMF_shannon, "Laurentian_Mixed_Forest"), FM(ME_LMF_shannon, "Laurentian_Mixed_Forest"), FM(MN_LMF_shannon, "Laurentian_Mixed_Forest"), FM(WI_LMF_shannon, "Laurentian_Mixed_Forest"), FM(MI_LMF_shannon, "Laurentian_Mixed_Forest"))

EF = bind_rows(MRF_shannon2, MF_shannon2, CMF_shannon2, EBF_shannon2, LMF_shannon2, EBFC_shannon2)
write.csv(EF, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\EF.csv")
EF$biome <- factor(EF$biome)
EF$lat <- as.numeric(EF$lat)

####Eastern Mountains

A_shannon2 <- bind_rows(FM(NY_A_shannon, "Adirondack"), FM(NE_A_shannon, "Adirondack"))
ABF_shannon2 <- bind_rows(FM(GA_ABF_shannon, "Appalachian Broadleaf Forest"), FM(TN_ABF_shannon, "Appalachian Broadleaf Forest"), FM(NC_ABF_shannon, "Appalachian Broadleaf Forest"), FM(VA_ABF_shannon, "Appalachian Broadleaf Forest"), FM(WV_ABF_shannon, "Appalachian Broadleaf Forest"), FM(MD_ABF_shannon, "Appalachian Broadleaf Forest"), FM(PA_ABF_shannon, "Appalachian Broadleaf Forest"),
                          FM(SC_ABF_shannon, "Appalachian Broadleaf Forest"))
EM = bind_rows(A_shannon2, ABF_shannon2)
write.csv(EM, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\EM.csv")
EM$biome <- factor(EM$biome)
EM$lat <- as.numeric(EM$lat)

#### TPP

MN_TPP_shannon <- gen_shannon(MN_TPP)
IA_TPP_shannon <- gen_shannon(IA_TPP)
MO_TPP_shannon <- gen_shannon(MO_TPP)
IL_TPP_shannon <- gen_shannon(IL_TPP)
NE_TPP_shannon <- gen_shannon(NE_TPP)
KS_TPP_shannon <- gen_shannon(KS_TPP)
OK_TPP_shannon <- gen_shannon(OK_TPP)

TPP_shannon2 <- bind_rows(FM(MN_TPP_shannon, "Temperate Prairie Parkland"), FM(IA_TPP_shannon, "Temperate Prairie Parkland"),FM(MO_TPP_shannon, "Temperate Prairie Parkland"), FM(IL_TPP_shannon, "Temperate Prairie Parkland"), FM(NE_TPP_shannon, "Temperate Prairie Parkland"), FM(KS_TPP_shannon, "Temperate Prairie Parkland"), FM(OK_TPP_shannon, "Temperate Prairie Parkland"))

#### GPS 
ND_GPS_shannon <- gen_shannon(ND_GPS)
SD_GPS_shannon <- gen_shannon(SD_GPS)
NE_GPS_shannon <- gen_shannon(NE_GPS)
KS_GPS_shannon <- gen_shannon(KS_GPS)
OK_GPS_shannon <- gen_shannon(OK_GPS)
CO_GPS_shannon <- gen_shannon(CO_GPS)
WY_GPS_shannon <- gen_shannon(WY_GPS)
MT_GPS_shannon <- gen_shannon(MT_GPS)

GPS_shannon2 <- bind_rows(FM(ND_GPS_shannon, "Great Plains Steppe"), 
                          FM(SD_GPS_shannon, "Great Plains Steppe"),
                          FM(NE_GPS_shannon, "Great Plains Steppe"),
                          FM(KS_GPS_shannon, "Great Plains Steppe"),
                          FM(OK_GPS_shannon, "Great Plains Steppe"),
                          FM(CO_GPS_shannon, "Great Plains Steppe"),
                          FM(WY_GPS_shannon, "Great Plains Steppe"),
                          FM(MT_GPS_shannon, "Great Plains Steppe"))
#### GPPS
WA_GPPS_shannon <- gen_shannon(WA_GPPS)
OR_GPPS_shannon <- gen_shannon(OR_GPPS)
NV_GPPS_shannon <- gen_shannon(NV_GPPS)
ID_GPPS_shannon <- gen_shannon(ID_GPPS)
WY_GPPS_shannon <- gen_shannon(WY_GPPS8)

GPPS_shannon2 <- bind_rows(FM(WA_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                           FM(OR_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                           FM(NV_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                           FM(ID_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                           FM(WY_GPPS_shannon, "Great Plains Palousse Dry Steppe"))
                           

#### SPP
TX_SPP_shannon <- gen_shannon(TX_SPP)
OK_SPP_shannon <- gen_shannon(OK_SPP)

SPP_shannon2 <- bind_rows(FM(OK_SPP_shannon, "Subtropical Prairie Parkland"), FM(TX_SPP_shannon, "Subtropical Prairie Parkland"))

G <- bind_rows(SPP_shannon2, GPPS_shannon2, GPS_shannon2, TPP_shannon2)
write.csv(G, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\G.csv")
G$biome <- factor(G$biome)
G$lat <- as.numeric(G$lat)

### Desert

CA_ASD_shannon <- gen_shannon(CA_ASD)
NV_ASD_shannon <- gen_shannon(NV_ASD)
AZ_ASD_shannon <- gen_shannon(AZ_ASD)

ASD_shannon2 <- bind_rows(FM(CA_ASD_shannon, "American Semi-Desert"), FM(NV_ASD_shannon, "American Semi-Desert"), FM(AZ_ASD_shannon, "American Semi-Desert"))

AZ_ADM_shannon <- gen_shannon(AZ_ADM)
NM_ADM_shannon <- gen_shannon(NM_ADM)

ADM_shannon2 <- bind_rows(FM(AZ_ADM_shannon, "Arizona-New Mexico Mountains Semi-Desert"), FM(NM_ADM_shannon, "Arizona-New Mexico Mountains Semi-Desert"))

AZ_CPAD_shannon <- gen_shannon(AZ_CPAD)
NM_CPAD_shannon <- gen_shannon(NM_CPAD7)
UT_CPAD_shannon <- gen_shannon(UT_CPAD)
CO_CPAD_shannon <- gen_shannon(CO_CPAD)

CPAD_shannon2 <- bind_rows(FM(AZ_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),FM(NM_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),FM(UT_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),FM(CO_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"))

AZ_CSD_shannon <- gen_shannon(AZ_CSD)
NM_CSD_shannon <- gen_shannon(NM_CSD)
TX_CSD_shannon <- gen_shannon(TX_CSD)

CSD_shannon2 <- bind_rows(FM(AZ_CSD_shannon, "Chihuahuan Semi-Desert"), FM(NM_CSD_shannon, "Chihuahuan Semi-Desert"), FM(TX_CSD_shannon, "Chihuahuan Semi-Desert"))

NM_SP_shannon <- gen_shannon(NM_SP)
TX_SP_shannon <- gen_shannon(TX_SP)

SP_shannon2 <- bind_rows(FM(NM_SP_shannon, "Southwest Plateau"),FM(TX_SP_shannon, "Southwest Plateau"))

NV_ISD_shannon <- gen_shannon(NV_ISD)
UT_ISD_shannon <- gen_shannon(UT_ISD)
CO_ISD_shannon <- gen_shannon(CO_ISD)

ISD_shannon2 <- bind_rows(FM(NV_ISD_shannon, "Inter-Mountain Semi-Desert"), FM(UT_ISD_shannon, "Inter-Mountain Semi-Desert"), FM(CO_ISD_shannon, "Inter-Mountain Semi-Desert"))


D <- bind_rows(ISD_shannon2, SP_shannon2, CSD_shannon2, CPAD_shannon2, ADM_shannon2, ASD_shannon2, ISD_shannon2)
write.csv(D, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\D.csv")
D$lat <- as.numeric(D$lat)

#### Western Forests

WA_PLMF_shannon <- gen_shannon(WA_PLMF)

OR_PLMF2_shannon <- gen_shannon(OR_PLMF2)
OR_PLMF4_shannon <- gen_shannon(OR_PLMF4)
OR_PLMF6_shannon <- gen_shannon(OR_PLMF6)
OR_PLMF10_shannon <- gen_shannon(OR_PLMF8)

OR_CCCF_shannon <- gen_shannon(OR_CCCF)
CA_CCCF_shannon <- gen_shannon(CA_CCCF)

CA_CCROW13_shannon <- gen_shannon(CA_CCROW13)
CA_CCROW3_shannon <- gen_shannon(CA_CCROW3)

WF <- bind_rows(FM(OR_CCCF_shannon, "California Coastal Chapparel Forest"),
                   FM(CA_CCCF_shannon, "California Coastal Chapparel Forest"),
                   FM(CA_CCROW13_shannon, "California Coastal Range Open Woodland"),
                   FM(CA_CCROW3_shannon, "California Coastal Range Open Woodland"))
                         
write.csv(WF, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\WF.csv")
WF$lat <- as.numeric(WF$lat)


#### Western Mountains
OR_CMF_shannon <- gen_shannon(OR_CMF)
WA_CMF_shannon <- gen_shannon(WA_CMF)
MT_NRMFS1_shannon <- gen_shannon(MT_NRMFS1)
MT_NRMFS2_shannon <- gen_shannon(MT_NRMFS2)
ID_NRMFS3_shannon <- gen_shannon(ID_NRMFS3)
ID_NRMFS4_shannon <- gen_shannon(ID_NRMFS4)
WA_NRMFS1_shannon <- gen_shannon(WA_NRMFS1)
WA_NRMFS2_shannon <- gen_shannon(WA_NRMFS2)
WA_NRMFS3_shannon <- gen_shannon(WA_NRMFS3)

WA_MRMS_shannon <- gen_shannon(WA_MRMS)
OR_MRMS4_shannon <- gen_shannon(OR_MRMS4)
OR_MRMS6_shannon <- gen_shannon(OR_MRMS6)
MT_MRMS_shannon <- gen_shannon(MT_MRMS)
ID_MRMS_shannon <- gen_shannon(ID_MRMS)

NV_NVUTM_shannon <- gen_shannon(NV_NVUTM)
UT_NVUTM6_shannon <- gen_shannon(UT_NVUTM6)
UT_NVUTM3_shannon <- gen_shannon(UT_NVUTM3)
CO_NVUTM_shannon <- gen_shannon(CO_NVUTM)

WY_SRMSOW1_shannon <- gen_shannon(WY_SRMSOW1)
WY_SRMSOW2_shannon <- gen_shannon(WY_SRMSOW2)
WY_SRMSOW4_shannon <- gen_shannon(WY_SRMSOW4)
CO_SRMSOW_shannon <- gen_shannon(CO_SRMSOW)
NM_SRMSOW8_shannon <- gen_shannon(NM_SRMSOW8)
NM_SRMSOW4_shannon <- gen_shannon(NM_SRMSOW4)
UT_SRMSOW_shannon <- gen_shannon(UT_SRMSOW)
ID_SRMSOW1_shannon <- gen_shannon(ID_SRMSOW1)


WM <- bind_rows(FM(OR_CMF_shannon, "Cascade Mixed Forest"),
                FM(WA_CMF_shannon, "Cascade Mixed Forest"),
                FM(MT_NRMFS1_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(MT_NRMFS2_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(ID_NRMFS3_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(ID_NRMFS4_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(WA_NRMFS1_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(WA_NRMFS2_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(WA_NRMFS3_shannon, "Northern Rocky Mountains Forest Steppe"),
                FM(WA_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                FM(OR_MRMS4_shannon, "Middle Rocky Mountain Steppe"),
                FM(OR_MRMS6_shannon, "Middle Rocky Mountain Steppe"),
                FM(MT_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                FM(ID_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                FM(NV_NVUTM_shannon, "Nevada-Utah Mountains"),
                FM(UT_NVUTM6_shannon, "Nevada-Utah Mountains"),
                FM(UT_NVUTM3_shannon, "Nevada-Utah Mountains"),
                FM(CO_NVUTM_shannon, "Nevada-Utah Mountains"),
                FM(WY_SRMSOW1_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(WY_SRMSOW2_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(WY_SRMSOW4_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(CO_SRMSOW_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(NM_SRMSOW8_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(NM_SRMSOW4_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(UT_SRMSOW_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                FM(ID_SRMSOW1_shannon, "Southern Rocky Mountain Steppe-Open Woodland"))

WM$lat <- as.numeric(WM$lat)


WM$lat <- as.numeric(WM$lat)
write.csv(WM, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\WM.csv")
all2 <- bind_rows(FM(OR_CCCF_shannon, "California Coastal Chapparel Forest"),
                  FM(CA_CCCF_shannon, "California Coastal Chapparel Forest"),
                  FM(CA_CCROW13_shannon, "California Coastal Range Open Woodland"),
                  FM(CA_CCROW3_shannon, "California Coastal Range Open Woodland"),
                  FM(OR_CMF_shannon, "Cascade Mixed Forest"),
                  FM(WA_CMF_shannon, "Cascade Mixed Forest"),
                  FM(MT_NRMFS1_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(MT_NRMFS2_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(ID_NRMFS3_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(ID_NRMFS4_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(WA_NRMFS1_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(WA_NRMFS2_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(WA_NRMFS3_shannon, "Northern Rocky Mountains Forest Steppe"),
                  FM(WA_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                  FM(OR_MRMS4_shannon, "Middle Rocky Mountain Steppe"),
                  FM(OR_MRMS6_shannon, "Middle Rocky Mountain Steppe"),
                  FM(MT_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                  FM(ID_MRMS_shannon, "Middle Rocky Mountain Steppe"),
                  FM(NV_NVUTM_shannon, "Nevada-Utah Mountains"),
                  FM(UT_NVUTM6_shannon, "Nevada-Utah Mountains"),
                  FM(UT_NVUTM3_shannon, "Nevada-Utah Mountains"),
                  FM(CO_NVUTM_shannon, "Nevada-Utah Mountains"), 
                  FM(NM_SP_shannon, "Southwest Plateau"),
                  FM(TX_SP_shannon, "Southwest Plateau"), 
                  FM(AZ_CSD_shannon, "Chihuahuan Semi-Desert"), 
                  FM(NM_CSD_shannon, "Chihuahuan Semi-Desert"), 
                  FM(TX_CSD_shannon, "Chihuahuan Semi-Desert"), 
                  FM(AZ_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),
                  FM(NM_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),
                  FM(UT_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"),
                  FM(CO_CPAD_shannon, "Colorado Plateau Semi-Arid Desert"), 
                  FM(AZ_ADM_shannon, "Arizona-New Mexico Mountains Semi-Desert"), 
                  FM(NM_ADM_shannon, "Arizona-New Mexico Mountains Semi-Desert"), 
                  FM(CA_ASD_shannon, "American Semi-Desert"), 
                  FM(NV_ASD_shannon, "American Semi-Desert"), 
                  FM(AZ_ASD_shannon, "American Semi-Desert"),
                  FM(ND_GPS_shannon, "Great Plains Steppe"), 
                  FM(SD_GPS_shannon, "Great Plains Steppe"),
                  FM(NE_GPS_shannon, "Great Plains Steppe"),
                  FM(KS_GPS_shannon, "Great Plains Steppe"),
                  FM(OK_GPS_shannon, "Great Plains Steppe"),
                  FM(CO_GPS_shannon, "Great Plains Steppe"),
                  FM(WY_GPS_shannon, "Great Plains Steppe"),
                  FM(MT_GPS_shannon, "Great Plains Steppe"),
                  FM(WA_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                  FM(OR_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                  FM(NV_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                  FM(ID_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                  FM(WY_GPPS_shannon, "Great Plains Palousse Dry Steppe"),
                  FM(OK_SPP_shannon, "Subtropical Prairie Parkland"), 
                  FM(TX_SPP_shannon, "Subtropical Prairie Parkland"),
                  FM(MN_TPP_shannon, "Temperate Prairie Parkland"), 
                  FM(IA_TPP_shannon, "Temperate Prairie Parkland"),
                  FM(MO_TPP_shannon, "Temperate Prairie Parkland"), 
                  FM(IL_TPP_shannon, "Temperate Prairie Parkland"), 
                  FM(NE_TPP_shannon, "Temperate Prairie Parkland"), 
                  FM(KS_TPP_shannon, "Temperate Prairie Parkland"), 
                  FM(OK_TPP_shannon, "Temperate Prairie Parkland"),
                  FM(PA_LMF_shannon, "Laurentian_Mixed_Forest"), 
                  FM(NY_LMF_shannon, "Laurentian_Mixed_Forest"), 
                  FM(ME_LMF_shannon, "Laurentian_Mixed_Forest"), 
                  FM(MN_LMF_shannon, "Laurentian_Mixed_Forest"), 
                  FM(WI_LMF_shannon, "Laurentian_Mixed_Forest"), 
                  FM(MI_LMF_shannon, "Laurentian_Mixed_Forest"),
                  FM(GA_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(TN_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(NC_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(VA_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(WV_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(MD_ABF_shannon, "Appalachian Broadleaf Forest"), 
                  FM(PA_ABF_shannon, "Appalachian Broadleaf Forest"),
                  FM(NY_A_shannon, "Adirondack"), 
                  FM(NE_A_shannon, "Adirondack"),
                  FM(PA_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(MA_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(TN_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(KY_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(WV_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(OH_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(NH_EBF_shannon, "Eastern Broadleaf Forest"), 
                  FM(ME_EBF_shannon, "Eastern Broadleaf Forest"),
                  FM(LA_MRF_shannon, "Mississippi River Forest"), 
                  FM(MI_MRF_shannon, "Mississippi River Forest"), 
                  FM(AK_MRF_shannon, "Mississippi River Forest"), 
                  FM(MO_MRF_shannon, "Mississippi River Forest"), 
                  FM(TN_MRF_shannon,"Mississippi River Forest"), 
                  FM(KY_MRF_Shannon, "Mississippi River Forest"),
                  FM(AL_MF_shannon, "Mixed Forest"), 
                  FM(AK_MF_shannon,  "Mixed Forest"), 
                  FM(LA_MF_shannon, "Mixed Forest"), 
                  FM(TX_MF_shannon,  "Mixed Forest"), 
                  FM(GA_MF_shannon,  "Mixed Forest"), 
                  FM(NC_MF_shannon, "Mixed Forest"), 
                  FM(VA_MF_shannon,  "Mixed Forest"), 
                  FM(SC_MF_shannon,  "Mixed Forest"), 
                  FM(MI_MF_shannon,  "Mixed Forest"),
                  FM(DA_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(MD_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(GA_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(VA_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(SC_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(NC_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(AL_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(MI_CMF_shannon, "Coastal Mixed Forest"), 
                  FM(LA_CMF_shannon, "Coastal Mixed Forest"),
                  FM(NV_ISD_shannon, "Inter-Mountain Semi-Desert"), 
                  FM(UT_ISD_shannon, "Inter-Mountain Semi-Desert"), 
                  FM(CO_ISD_shannon, "Inter-Mountain Semi-Desert"),
                  FM(EBFC_shannon, "Eastern Broadleaf Forest (Continental)"),
                  FM(WY_SRMSOW1_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(WY_SRMSOW2_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(WY_SRMSOW4_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(CO_SRMSOW_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(NM_SRMSOW8_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(NM_SRMSOW4_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(UT_SRMSOW_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(ID_SRMSOW1_shannon, "Southern Rocky Mountain Steppe-Open Woodland"),
                  FM(AR_SMF_shannon, "Southeastern Mixed Forest"),
                  FM(OK_SMF_shannon, "Southeastern Mixed Forest"))

write.csv(all2, "C:\\Users\\andre\\Documents\\BIOL7800\\Final Project\\Shannon Only\\all.csv")                  
all2$biome <- factor(all2$biome)
all2$biome <- relevel(all2$biome, ref = "Southeastern Mixed Forest")

all2$lat <- as.numeric(all2$lat)
all2

SLmodelWM <- gen_m0(WM)
modelWM <- aov(Shannon~lat+biome+lat*biome, data = WM)
ANOVA_WM<-anova(modelWM)
ANOVA_WM
Residuals_WM<- residuals(modelWM)
qqnorm(Residuals_WM, main="QQ Plot")
qqline(Residuals_WM, col = 2) 
modelWMLat <- gen_m1(WM)
coef(modelWMLat)[2]


SLmodelWF<-gen_m0(WF)
modelWF <- aov(Shannon~lat+biome+lat*biome, data = WF)
ANOVA_WF<-anova(modelWF)
ANOVA_WF
modelWFLat <- gen_m1(WF)
coef(modelWFLat)[2]



SLmodelD <- gen_m0(D)
modelD <- aov(Shannon~lat+biome+lat*biome, data = D)
ANOVA_D<-anova(modelD)
ANOVA_D
modelDLat <- gen_m1(D)
coef(modelDLat)[2]


SLmodelG<-gen_m0(G)
modelG <- aov(Shannon~lat+biome+lat*biome, data = G)
ANOVA_G<-anova(modelG)
ANOVA_G
modelGLat <- gen_m1(G)
coef(modelGLat)[2]

SLmodelEM<-gen_m0(EM)
modelEM <- aov(Shannon~lat+biome+lat*biome, data = EM)
ANOVA_EM<-anova(modelEM)
ANOVA_EM
Residuals_EM<- residuals(modelEM)
qqnorm(Residuals_EM, main="QQ Plot")
qqline(Residuals_EM, col = 2) 
modelEMLat <- gen_m1(EM)
coef(modelEMLat)[2]

SLmodelEF<- gen_m0(EF)
modelEF <- aov(Shannon~lat+biome+lat*biome, data = EF)
ANOVA_EF<-anova(modelEF)
ANOVA_EF
Residuals_EF <- residuals(modelEF)
qqnorm(Residuals_EF, main="QQ Plot")
qqline(Residuals_EF, col = 2) 
modelEFLat <- gen_m1(EF)
coef(modelEFLat)[2]


modelall <- gen_m0(all2)
model1 <- aov(Shannon~lat+biome+lat*biome, data = all2)
Residuals_all <- residuals(model1)
ANOVA_all <- anova(model1)
ANOVA_all
qqnorm(Residuals_all, main="QQ Plot All")
qqline(Residuals_all, col = 2) 
modelall
model <- aov(Shannon ~ lat+biome+lat*biome, data = all2)
plot(model, 1)
modelallLat <- gen_m1(all2)
coef(modelallLat)[2]
plot(model, 3)

plotWM <- ggplot(WM, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 5 Western Mountain Biomes")

plotWF <- ggplot(WF, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 3 Western Forest Biomes")

plotD <- ggplot(D, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 7 Desert Biomes")

plotEM <- ggplot(EM, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 2 Eastern Mountain Biomes")

plotEF <- ggplot(EF, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 6 Eastern Forest Biomes")

plotG <- ggplot(G, aes(x= lat, y = Shannon, color = biome)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, aes(group=biome)) +
  xlab("Latitude") + 
  ylab("Shannon Diversity Index") +
  ggtitle("The Impact of Latitude on the Shannon Diversity Index Values of 4 Grassland Biomes")

### Jaccard

# Create a list of 24 data frames (replace this with your actual list of data frames)

EBF <- mutate(EBF, biome = rep(c("Eastern Broadleaf Forest")))
Mississippi_River_Forest <- mutate(Mississippi_River_Forest, biome = rep(c("Mississippi River Forest")))
Mixed_Forest <- mutate(Mixed_Forest, biome = rep(c("Mixed Forest")))
Coastal_Mixed_Forest <- mutate(Coastal_Mixed_Forest, biome = rep(c("Coastal Mixed Forest")))
LMF <- mutate(LMF, biome = rep(c("Laurentian Mixed Forest")))
SMF <- mutate(SMF, biome = rep(c("Southeastern Mixed Forest")))

write.csv(CCCF, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CCCF.csv")
write.csv(CA_CCROW, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CA_CCROW.csv")
write.csv(NRMFS, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\NRMFS.csv")
write.csv(CAMF, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CAMF.csv")
write.csv(MRMS, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\MRMS.csv")
write.csv(NVUTM, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\NVUTM.csv")
write.csv(SP, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\SP.csv")
write.csv(CSD, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CSD.csv")
write.csv(CPAD, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CPAD.csv")
write.csv(ADM, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\ADM.csv")
write.csv(ASD, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\ASD.csv")
write.csv(GPS, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\GPS.csv")
write.csv(GPPS, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\GPPS.csv")
write.csv(SPP, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\SPP.csv")
write.csv(TPP, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\TPP.csv")

write.csv(Appalacian_Broadleaf_Forest, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\ABF.csv")
write.csv(Adirondack, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\Adirondack.csv")
write.csv(ISD, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\ISD.csv")
write.csv(EBFC, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBFC.csv")
write.csv(SRMSOW, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\SRMSOW.csv")

EBFC <- mutate(EBFC, biome = rep(c("Eastern Broadleaf Forest (Continental)")))
SRMSOW <- mutate(SRMSOW, biome = rep(c("Southern Rocky Mountain Steppe-Open Woodland")))

split_size <- nrow(LMF) / 2
LMF1 <- LMF[1:split_size, ]
LMF2 <- LMF[(split_size + 1):nrow(LMF), ]
write.csv(LMF, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\LMF1.csv")
write.csv(LMF, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\LMF2.csv")

split_size <- nrow(EBFC) / 3
EBFC1 <- EBFC[1:split_size, ]
EBFC2 <- EBFC[(split_size + 1):(2 * split_size), ]
EBFC3 <- EBFC[(2 * split_size + 1):nrow(EBFC), ]
write.csv(EBFC1, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBFC1.csv")
write.csv(EBFC2, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBFC2.csv")
write.csv(EBFC3, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBFC3.csv")

split_size <- nrow(Mixed_Forest) / 3
MF1 <- Mixed_Forest[1:split_size, ]
MF2 <- Mixed_Forest[(split_size + 1):(2 * split_size), ]
MF3 <- Mixed_Forest[(2 * split_size + 1):nrow(Mixed_Forest), ]
write.csv(MF1, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\MF1.csv")
write.csv(MF2, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\MF2.csv")
write.csv(MF3, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\MF3.csv")

split_size <- nrow(Coastal_Mixed_Forest) / 3
CMF1 <- Coastal_Mixed_Forest[1:split_size, ]
CMF2 <- Coastal_Mixed_Forest[(split_size + 1):(2 * split_size), ]
CMF3 <- Coastal_Mixed_Forest[(2 * split_size + 1):nrow(Coastal_Mixed_Forest), ]
write.csv(CMF1, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CMF1.csv")
write.csv(CMF2, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CMF2.csv")
write.csv(CMF3, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\CMF3.csv")

split_size <- nrow(EBF) / 4
EBF1 <- EBF[1:split_size, ]
EBF2 <- EBF[(split_size + 1):(2 * split_size), ]
EBF3<- EBF[(2 * split_size + 1):(3 * split_size), ]
EBF4 <- EBF[(3 * split_size + 1):nrow(EBF), ]
write.csv(EBF1, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBF1.csv")
write.csv(EBF2, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBF2.csv")
write.csv(EBF3, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBF3.csv")
write.csv(EBF4, "C:\\Users\\andre\\Documents\\BIOL7800\\Project Proposal\\Cleaned_Data\\EBF4.csv")
