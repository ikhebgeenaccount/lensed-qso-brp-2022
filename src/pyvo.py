import pyvo
url = "https://irsa.ipac.caltech.edu/SCS?table=fp_psc"

ra = 240.41854
dec = 43.27994

print(pyvo.conesearch(url, pos=(ra, dec), radius=.09))
