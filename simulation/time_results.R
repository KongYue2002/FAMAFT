### computation time ###
time_famaft <- mean(time.famaft.N)
time_maft <- mean(time.maft.N)
time_tf <- mean(time.tf.N)
time <- data.frame(FAMAFT = time_famaft, MAFT = time_maft, TF = time_tf) 

write.csv(time, file = "time.csv", row.names = TRUE)