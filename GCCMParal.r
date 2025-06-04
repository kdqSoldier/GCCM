#  处理每一个窗口，计算交叉映射
GCCMSingle<-function(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,cores=NULL)
{
  print("GCCMSingle")
  # 创建一个Data frame
  x_xmap_y <- data.frame()
  # 如果没有指定核心数，单线程执行
  if(is.null(cores))
  {
    # 遍历每一行
    for(r in 1:(totalRow-lib_size+1))
    {
      # 调用 GCCMSingleInner 进行具体的计算
      x_xmap_y<-rbind(x_xmap_y,GCCMSingleInner(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r))
    }
  # 如果指定了核心数，多线程并行执行
  }else
  {
    # 创建并行集群
    cl <- makeCluster(cores)
    # 创建并行集群
    registerDoParallel(cl)
    # 创建并行集群
    clusterExport(cl,deparse(substitute(GCCMSingleInner)))
    clusterExport(cl,deparse(substitute(locate)))
    clusterExport(cl,deparse(substitute(projection)))
    clusterExport(cl,deparse(substitute(distance_Com)))
    clusterExport(cl,deparse(substitute(compute_stats)))
    #browser()
    # 这里
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print(lib_size)
    print(b)
    print(cores)
    # 创建并行集群
    x_xmap_y <- foreach(r=(1:(totalRow-lib_size+1)), .combine='rbind') %dopar% GCCMSingleInner(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r)
    # print(x_xmap_y)
    # browser()
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 停止并行集群
    stopCluster(cl)
  }

  # browser()
  return(x_xmap_y)
}

# 处理每一个嵌入维度的具体计算
GCCMSingleInner<-function(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,r)
{

  browser()
  x_xmap_y <- data.frame()
  lib_ids <- NULL
  #if (lib_size==80) {
  #  x_xmap_y <- 99999999
  #}
  # 遍历列
  for(c in 1:(totalCol-lib_size+1))
  {
    
    # browser()
    # 初始化预测索引
    pred_indices <- rep.int(FALSE, times = totalRow*totalCol)
    # 初始化库索引
    lib_indices<- rep.int(FALSE, times = totalRow*totalCol)
    # 设置预测索引为 TRUE
    pred_indices[locate(pred[,1],pred[,2],totalRow,totalCol) ]<-TRUE
    # 将 NA 值对应的索引设置为 FALSE
    pred_indices[which(is.na(yPred)) ]<-FALSE
    
    # 库行序列
    lib_rows<-seq(r,(r+lib_size-1))
    # 库列序列
    lib_cols<-seq(c,(c+lib_size-1))
    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # 合并库行列
    lib_ids<-merge(lib_rows,lib_cols)
    #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # 设置库索引为 TRUE
    lib_indices[locate(lib_ids[,1],lib_ids[,2],totalRow,totalCol)]<-TRUE
    
    # 如果目标矩阵的线性窗口中的缺失值超过一半，跳过
    if(length(which(is.na(yPred[which(lib_indices)]))) > ((lib_size*lib_size)/2))
    {
      #x_xmap_y <- 99999999999
      next
    }
    # 运行交叉映射并存储结果
    # run cross map and store results
    results <-  projection(xEmbedings,yPred,lib_indices ,pred_indices,b)
    # 将结果添加到数据框中
    x_xmap_y <- rbind(x_xmap_y, data.frame(L = lib_size, rho = results$stats$rho)) 
  }
  
  return(x_xmap_y)
  # return(lib_ids)
}


# 主函数，处理整个图像矩阵并调用 GCCMSingle 进行计算
GCCM<-function(xMatrix, yMatrix, lib_sizes, lib, pred, E, b, tau = 1, cores=NULL)
{

  print("GCCM")
  #print(tau)
  #print(b)
  #print(cores)
  # 获取图像矩阵的尺寸
  imageSize<-dim(xMatrix)
  # 总行数
  totalRow<-imageSize[1]
  # 总列数
  totalCol<-imageSize[2]
  
  # 转置 y 矩阵并转换为数组
  yPred<- as.array(t(yMatrix))
  
  xEmbedings<-list()
  # 转置 x 矩阵并转换为数组
  xEmbedings[[1]]<- as.array(t(xMatrix))
  
  for(i in 1:E)
  {
    # 构建滞后变量的嵌入矩阵
    xEmbedings[[i+1]]<-laggedVariableAs2Dim(xMatrix, i)  #### row first
  
  }
  
  x_xmap_y <- data.frame()
  # 遍历每一个库大小
  for(lib_size in lib_sizes)
  {
    # 调用 GCCMSingle 进行计算，并合并结果
     x_xmap_y<-rbind(x_xmap_y,GCCMSingle(xEmbedings,yPred,lib_size,pred,totalRow,totalCol,b,cores))

  }
  
  return (x_xmap_y)
}

# 将二维坐标转换为一维索引
locate<-function(curRow,curCOl,totalRow,totalCol)
{
  # 遍历每一个库大小
  return ((curRow-1)*totalCol+curCOl) 
}
# 进行预测，计算权重并进行加权预测
projection<-function(embedings,target,lib_indices, pred_indices,num_neighbors)
{
  print("projection")
  # 初始化预测结果为 NaN
  pred <- rep.int(NaN, times = length(target))
  # 遍历每一个预测索引
  for(p in which (pred_indices))
  { 
    # 临时存储当前库索引状态
    temp_lib <- lib_indices[p]
    # 将当前库索引设置为 FALSE
    lib_indices[p] <- FALSE
    # 获取所有库索引
    libs <- which(lib_indices)
    # 计算嵌入向量之间的距离
    distances<-distance_Com(embedings,libs,p)
    
    #distances<-colMeans(distances)
    
    # find nearest neighbors
    # 找到最近邻
    neighbors <- order(distances)[1:num_neighbors]
    min_distance <- distances[neighbors[1]]
    if(is.na(min_distance))
    {
      # 如果最小距离为 NaN，跳过
      lib_indices[p] <- temp_lib 
      next
    }
    # compute weights
    # 计算权重
    if(min_distance == 0) # perfect match
    {
      weights <- rep.int(0.000001, times = num_neighbors)
      weights[distances[neighbors] == 0] <- 1
    }
    else
    {
      weights <- exp(-distances[neighbors]/min_distance)
      weights[weights < 0.000001] <- 0.000001
    }
    total_weight <- sum(weights)
    
    # make prediction
    # 进行预测
    pred[p] <- (weights %*% target[libs[neighbors]]) / total_weight
    
    # 恢复库索引状态
    lib_indices[p] <- temp_lib 
  }
  
  # return output & stats
  # 返回结果和统计信息
  return(list(pred = pred, stats = compute_stats(target[pred_indices], pred[pred_indices])))
  
}

# 计算嵌入向量之间的距离
distance_Com<-function(embeddings,libs,p)
{
  print("distance_Com")
  # 初始化距离向量
  distances<-c()
  # 遍历每一个嵌入向量
  for(e in 1:length(embeddings))
  {
    emd<-embeddings[[e]]
    # 复制当前点的嵌入向量
    q <- matrix(rep(emd[p], length(libs)), nrow = length(libs), byrow = T)
  
    # 计算绝对距离并合并到距离向量中
    distances<-cbind(distances,abs(emd[libs]-emd[p]))
    
  }
  # 返回距离向量的行均值
  return (rowMeans(distances,na.rm=TRUE))
  
}

