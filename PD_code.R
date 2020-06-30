y = as.factor(y)

pd_check = check_data(x = data, y = y, B = 100, k = 10 , repB = 2, methods = c("boruta","RFE"), topN = 20)
