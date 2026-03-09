# cube root function (~ log-transformation with negative data)       
nthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}