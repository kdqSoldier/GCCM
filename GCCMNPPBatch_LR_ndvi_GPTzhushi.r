# 设置工作目录
setwd("D:/ProJect/Rstudio/NPP/1")

# 加载并行计算所需的包
library(parallel)
library(foreach)
library(doParallel)
library(iterators)

# 加载自定义函数脚本（需包含 GCCM、significance、confidence 等函数）
source("basic.r")
source("GCCMParal.r")

# 加载 NDVI 数据和气候变量图像（climateImages）
load("ndvi_20_10000_2.RData")

# ----------------------------------------------------
# 变量初始化
# ----------------------------------------------------
# 设置候选解释变量名称列表（这里使用20年尺度的变量）
xNames <- c("DEM_20_10000", "dh_01_20_10000", "landuse_01_20_10000", "prc_01_20_10000",
            "rh_01_20_10000", "sheep_01_20_10000", "slope_20_10000", 
            "soil_20_10000", "tem_01_20_10000", "vege_20_10000")

# yName：响应变量的图层名称（NDVI）
yName <- "ndvi"

# 嵌入维度（E）：用于构造状态空间的维度(CCM方法中常见取值为2~4)
E <- 4
# 嵌入步长：延迟嵌入时间间隔(通常为 E+1)         
b <- E + 2  
# lib_sizes：交叉映射分析窗口大小序列(通常设置为 10~最大观测点数量的均匀序列)
lib_sizes <- seq(40, 65, 5)

# 设置预测点过滤参数,稀疏采样使用（用于降低计算负担）
# 起始像元索引（行列均为5）
boundary_XY <- 5
# 采样步长（跳2个像素选一个点）
boundary_BY <- 2

# ----------------------------------------------------
# 构建预测网格与空间坐标
# ----------------------------------------------------
# 从NDVI图像中提取维度信息（行数、列数）
yMatrix <- ndviImage
# imageSize = [行数，列数]
imageSize <- dim(yMatrix)
totalRow <- imageSize[1]
totalCol <- imageSize[2]

# 构建稀疏预测坐标点网格，稀疏采样（加速计算）
predRows <- seq(boundary_XY, totalRow, boundary_BY)
predCols <- seq(boundary_XY, totalCol, boundary_BY)
# 生成所有预测点组合（row, col）
pred <- merge(predRows, predCols)

# 构建整个图像的像元坐标矩阵 coords，用于趋势建模
coodsX <- seq(1, totalRow)
coodsY <- seq(1, totalCol)
# 形成 full grid 坐标表
coords <- merge(coodsX, coodsY)
colnames(coords) <- c("coordX", "coodY")

# ----------------------------------------------------
# 响应变量 NDVI 去趋势
# ----------------------------------------------------
# 对 NDVI 图像数据去趋势处理（消除空间趋势）
# 将 NDVI 图像展开为一维向量
y <- as.vector(yMatrix)
# 使用线性模型 y ~ x + y 坐标，拟合并剔除空间趋势项
lmModel <- lm(y ~ coordX + coodY, data = cbind(y, coords))
# 计算去趋势后的 NDVI 值（残差）
y_detrended <- y - predict(lmModel, coords)
# 重新转为矩阵形式，得到去趋势后的响应变量图像
yMatrixM <- matrix(y_detrended, nrow = totalRow, ncol = totalCol)

# ----------------------------------------------------
# 遍历候选解释变量进行 CCM 分析
# ----------------------------------------------------
# 此处二选一即可：：
# 选取参数分析，只计算第4个变量
for (c in seq(3, 3)) 
# 遍历每一个候选解释变量
# for(c in seq(1:length(climateImages)))
{
  
  # 提取当前变量图像
  xMatrix <- climateImages[[c]]
  
  # 对该解释变量图像也做空间趋势去除，和65行相同方法
  x <- as.vector(xMatrix)
  lmModel <- lm(x ~ coordX + coodY, data = cbind(x, coords))
  x_detrended <- x - predict(lmModel, coords)
  xMatrixM <- matrix(x_detrended, nrow = totalRow, ncol = totalCol)
  
  # 记录计算起始时间
  startTime <- Sys.time()
  
  # 执行 Convergent Cross Mapping 分析
  # x_xmap_y：表示 x预测y 的结果
  # y_xmap_x：表示 y预测x 的结果（用于判断因果方向）
  x_xmap_y <- GCCM(xMatrixM, yMatrixM, lib_sizes, lib = NULL, pred, E, b, cores = 32)
  y_xmap_x <- GCCM(yMatrixM, xMatrixM, lib_sizes, lib = NULL, pred, E, b, cores = 32)
  
  # 记录计算结束时间并输出运行时间
  endTime <- Sys.time()
  print(difftime(endTime, startTime, units = "mins"))
  
  # 计算不同 library size 下的 rho（预测精度）平均值
  x_xmap_y$L <- as.factor(x_xmap_y$L)
  x_xmap_y_means <- do.call(rbind, lapply(split(x_xmap_y, x_xmap_y$L), function(x) {
    max(0, mean(x$rho, na.rm = TRUE))
  }))
  
  y_xmap_x$L <- as.factor(y_xmap_x$L)
  y_xmap_x_means <- do.call(rbind, lapply(split(y_xmap_x, y_xmap_x$L), function(x) {
    max(0, mean(x$rho, na.rm = TRUE))
  }))
  
  # 确保 lib_sizes 长度一致，避免维度不匹配
  lib_sizes <- lib_sizes[1:length(x_xmap_y_means)]
  
  # 预测点索引定位，预测点位置索引（用于提取观测值）
  predIndices <- locate(pred[,1], pred[,2], totalRow, totalCol)
  yPred <- as.array(t(yMatrix))
  predicted <- na.omit(yPred[predIndices])
  
  
  # ----------------------------------------------------
  # 显著性与置信区间计算
  # ----------------------------------------------------
  # 计算显著性和置信区间，检验预测显著性（与随机情况对比）
  x_xmap_y_Sig <- significance(x_xmap_y_means, length(predicted))
  y_xmap_x_Sig <- significance(y_xmap_x_means, length(predicted))
 
  # 计算置信区间（如95%）
  x_xmap_y_interval <- confidence(x_xmap_y_means, length(predicted))
  colnames(x_xmap_y_interval) <- c("x_xmap_y_upper", "x_xmap_y_lower")
  
  y_xmap_x_interval <- confidence(y_xmap_x_means, length(predicted))
  colnames(y_xmap_x_interval) <- c("y_xmap_x_upper", "y_xmap_x_lower")
  
  # ----------------------------------------------------
  # 保存结果到 CSV 和图像
  # ----------------------------------------------------
  # 汇总结果为数据框
  results <- data.frame(lib_sizes,x_xmap_y_means, y_xmap_x_means, x_xmap_y_Sig, y_xmap_x_Sig, x_xmap_y_interval, y_xmap_x_interval)
  
  # 输出中间结果（可选）
  print(x_xmap_y_means)
  print(y_xmap_x_means)
  
  # 保存为 CSV 文件
  xName <- xNames[c]
  filecsv <- paste("result_CSV_", xName, "_", yName, "_Pral.csv", sep = "")
  print(filecsv)
  write.csv(results, file = filecsv)
  
  # 绘图保存为 JPEG 图像
  filejpg <- paste("result_JPG_", xName, "_", yName, "_Pral.jpg", sep = "")
  print(filejpg)
  
  # 绘图展示 rho 随 lib_size 的变化趋势（验证因果方向）
  jpeg(filename = filejpg, width = 1300, height = 800)
  plot(lib_sizes, x_xmap_y_means, type = "l", col = "royalblue", lwd = 2,
       xlim = c(min(lib_sizes), max(lib_sizes)), ylim = c(0.0, 1),
       xlab = "L", ylab = "", cex.lab = 2, cex.axis = 2)
  mtext(expression(rho), side = 2, line = 2.5, cex = 2.3)
  lines(lib_sizes, y_xmap_x_means, col = "red3", lwd = 2)
  legend(min(lib_sizes), 1, legend = c("x xmap y", "y xmap x"),
         xjust = 0, yjust = 1, lty = 1, lwd = 2,
         col = c("royalblue", "red3"), cex = 2)
  dev.off()
  
  # 控制台输出运行信息
  cat("***------------------------***\n")
  cat("数据序号: ", c, "\n")
  cat("E: ", E, "\n")
  cat("b: ", b, "\n")
  cat("boundary_XY: ", boundary_XY, "\n")
  cat("boundary_BY: ", boundary_BY, "\n")
  cat("***------------------------***\n")
}
