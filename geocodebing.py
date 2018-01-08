import geocoder
import csv
import os

directory = 'ARB/'
outdirectory = 'ARB_OUT/'
for i in range(1, 59):

    with open(directory + str(i) + ".csv", 'rb') as csvfile, open(outdirectory + str(i) + '.csv', 'a') as csvout:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames + ['lat'] + ['lon']  # add new columns
        writer = csv.DictWriter(csvout, fieldnames)
        writer.writeheader()

        for row in reader:
            address = row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP']
            g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
            newrow = dict(row)
            if g.latlng: # g.json.has_key('lat'):
                # print("1")
                # print(g.json['lat'], g.json['lng'])
                newrow['lat'] = g.json['lat']
                newrow['lon'] = g.json['lng']
                writer.writerow(newrow) # Only write row if successfully geocoded
            else:
                address = row['FSTREET'] + ', ' + row['FCITY'] + ', California'
                g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                newrow = dict(row)
                if g.latlng:
                    print("2")
                    print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                    print(g.json['lat'], g.json['lng'])
                    newrow['lat'] = g.json['lat']
                    newrow['lon'] = g.json['lng']
                    writer.writerow(newrow) # Only write row if successfully geocoded
                else:
                    address = row['FSTREET'] + ', California'
                    g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                    newrow = dict(row)
                    if g.latlng:
                        print("3")
                        print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                        print(g.json['lat'], g.json['lng'])
                        newrow['lat'] = g.json['lat']
                        newrow['lon'] = g.json['lng']
                        writer.writerow(newrow) # Only write row if successfully geocoded
                    else:
                        address = row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP']
                        g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                        if g.latlng:
                            print("4")
                            print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                            print(g.json['lat'], g.json['lng'])
                            newrow['lat'] = g.json['lat']
                            newrow['lon'] = g.json['lng']
                            writer.writerow(newrow) # Only write row if successfully geocoded
                        else:
                            address = row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California'
                            g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                            if g.latlng:
                                print("5")
                                print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                                print(g.json['lat'], g.json['lng'])
                                newrow['lat'] = g.json['lat']
                                newrow['lon'] = g.json['lng']
                                writer.writerow(newrow) # Only write row if successfully geocoded
                            else:
                                address = row['FNAME'] + ', ' + row['FSTREET'] + ', California'
                                g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                                if g.latlng:
                                    print("6")
                                    print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                                    print(g.json['lat'], g.json['lng'])
                                    newrow['lat'] = g.json['lat']
                                    newrow['lon'] = g.json['lng']
                                    writer.writerow(newrow) # Only write row if successfully geocoded
                                else:
                                    address = row['FNAME'] + ', California'
                                    g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                                    if g.latlng:
                                        print("7")
                                        print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                                        print(g.json['lat'], g.json['lng'])
                                        newrow['lat'] = g.json['lat']
                                        newrow['lon'] = g.json['lng']
                                        writer.writerow(newrow) # Only write row if successfully geocoded
                                    else:
                                        address = row['FNAME']
                                        g = geocoder.bing(address, key='xpBPjeRgvbu9W7K7ZDkW~xf_tJOQYrZJJa6kc8DQHaQ~AnLiS2LtDjp1w2CKrnQ_evZj0QFI0Qot6nXwiTa_VIp7E8JyP-ztrfgk6Q9Qk7-k')
                                        if g.latlng:
                                            print("8")
                                            print(row['FNAME'] + ', ' + row['FSTREET'] + ', ' + row['FCITY'] + ', California ' + row['FZIP'])
                                            print(g.json['lat'], g.json['lng'])
                                            newrow['lat'] = g.json['lat']
                                            newrow['lon'] = g.json['lng']
                                            writer.writerow(newrow) # Only write row if successfully geocoded
