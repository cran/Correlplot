gfunc <-
function (z, x, u, v, w) {
  return ((z + (w / outer (u, v)) * (x - z)) * sqrt(outer(u, v)))
}
