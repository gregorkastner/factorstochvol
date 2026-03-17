# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
  .Call(`_factorstochvol_RcppExport_registerCCallable`)
})