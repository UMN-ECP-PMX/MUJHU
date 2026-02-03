
anyNA(out)

lt0 <- function(var){
  var <- sym(var)
  out %>% count(!!var < 0) 
}

lt0("depot")
lt0("cent")
lt0("periph")
lt0("centv")
lt0("cente")
lt0("centr")
lt0("centvdp")
lt0("centedp")
lt0("centrdp")

lt0("trans1")
lt0("trans2")
lt0("trans3")
lt0("trans4")
lt0("trans5")
lt0("trans6")
lt0("trans7")

out %>% filter(depot < 0)

