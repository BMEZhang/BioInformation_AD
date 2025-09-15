vis_GOKEGG_zwr <- function(val,ego_bp,kk){
  result = switch(  
    val,  
    "barplot"= {print("test")
      print("barplot created")},  
    "dotplot"= {print("test")
      print("dotplot created")},  
    "ego_bp2"= {print("test")
      print("ego_bp2 created")},  
    "cnetplot"= {print("test")
      print("cnetplot created")}, # 显示所有节点标签
    print('no graph')
  )
}


vis_opt = 2
vis_GOKEGG_zwr(vis_opt,vis_opt,vis_opt)

