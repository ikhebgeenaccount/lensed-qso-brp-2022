import pyvo
url = "https://irsa.ipac.caltech.edu/SCS?table=fp_psc"

ra = 240.41854
dec = 43.27994

print(pyvo.conesearch(url, pos=(ra, dec), radius=.09))


# Spitzer stuff
# VO simple conesearch
# https://irsa.ipac.caltech.edu/docs/vo_scs.html
# Catalog search api
# https://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html
# Search query result
# https://irsa.ipac.caltech.edu/workspace/TMP_MMhDTh_25060/Gator/irsa/20890/tbview.html
# Gator catalog list
# https://irsa.ipac.caltech.edu/applications/Gator/
