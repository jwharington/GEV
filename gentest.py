from scipy.stats import genextreme

# genextreme.pdf(xi, mu, sigma)
r = genextreme.rvs(0.5, 0.2, 0.3, size=5000)
for it in r:
    print(it)

# returns: mu, sigma, xi
#
#     .1991     .2954     .5047   CONVGD
