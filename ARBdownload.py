from urllib import request
import pandas as pd

# download all 58 counties' data
def download(outdirectory):
    for i in range(1, 59):
    # for i in range(19, 20):
        request.urlretrieve("https://www.arb.ca.gov/app/emsinv/facinfo/faccrit_output.csv?&dbyr=2015&ab_=&dis_=&co_=" + str(i) + "&fname_=&city_=&sort=FacilityNameA&fzip_=&fsic_=&facid_=&all_fac=C&chapis_only=&CERR=&dd=", outdirectory + str(i) + ".csv")

# download the extra data for each factility in a county
def download_facility(infile, outdirectory):
    facilities = pd.read_csv(infile)
    for index, row in facilities.iterrows():
        # temp = request.urlretrieve("https://www.arb.ca.gov/app/emsinv/facinfo/facdet_output.csv?&dbyr=2015&ab_=MC&dis_=AMA&co_=3&fname_=&city_=&sort=C&fzip_=&fsic_=&facid_=" + str(row['FACID']) + "&all_fac=&chapis_only=&CERR=&dd=")
        temp_f = pd.read_csv("https://www.arb.ca.gov/app/emsinv/facinfo/facdet_output.csv?&dbyr=2015&ab_=" + str(row['AB']) + "&dis_=" + str(row['DIS']) + "&co_=" + str(row['CO']) + "&fname_=&city_=&sort=C&fzip_=&fsic_=&facid_=" + str(row['FACID']) + "&all_fac=&chapis_only=&CERR=&dd=")
        facilities.set_value(index, 'PM2.5T', temp_f.loc[temp_f['POLLUTANT NAME'] == 'PM2.5'].iloc[0]['EMISSIONS_TONS_YR'])
    facilities.to_csv('test.csv')


def main():
    outdirectory = 'ARB/'
    # download(outdirectory)
    outdirectory = 'fresno_facilities'
    download_facility('fresno_industry_arb.csv', outdirectory)




if __name__ == "__main__":
    main()
