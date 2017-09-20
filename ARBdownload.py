import urllib

outdirectory = 'ARB/'

# go through all 58 counties
for i in range(1, 59):
# for i in range(1, 2):
    urllib.urlretrieve("https://www.arb.ca.gov/app/emsinv/facinfo/faccrit_output.csv?&dbyr=2015&ab_=&dis_=&co_=" + str(i) + "&fname_=&city_=&sort=FacilityNameA&fzip_=&fsic_=&facid_=&all_fac=C&chapis_only=&CERR=&dd=", outdirectory + str(i) + ".csv")
