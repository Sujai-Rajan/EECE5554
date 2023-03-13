#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import utm
import serial 
import rospy
from gnss_driver.msg import *
from gnss_driver.msg import gnss_msg
from datetime import datetime

def latlon_conv(loc):
    deg = int(loc/100)
    conv = deg+(loc/100-deg)*10/6
    return conv

def gnss_driver():
    rospy.init_node('gnss',anonymous=True)
    gnss_data = rospy.Publisher('/gnss', gnss_msg, queue_size=10)
    rate = rospy.Rate(10) #10Hz

    connected_port = sys.argv[1]
    serial_port = rospy.get_param('~port',connected_port)
    serial_baud = rospy.get_param('~baudrate',4800)

    sport = serial.Serial(serial_port,serial_baud,timeout = 3)
    print(serial_port)
    print('\n')

    while not rospy.is_shutdown():
        full_gnss_str = sport.readline()
        #full_gnss_str = full_gnss_str.decode('utf-8')

        if full_gnss_str == '':
            rospy.logwarn("No data from RTK GNSS.")
            print("No data from RTK GNSS")
        else:
            if full_gnss_str.startswith(b"\r$GNGGA") or full_gnss_str.startswith(b"$GNGGA"):
                full_gnss_str = full_gnss_str.decode()
                gnss_str = full_gnss_str.split(",")
                print(gnss_str)
                print('\n')

                utc_str = gnss_str[1]
                utc_timestamp = datetime.strptime(utc_str,"%H%M%S.%f")
                utc_hrs = utc_timestamp.hour
                utc_min = utc_timestamp.minute
                utc_sec = utc_timestamp.second
                utc_msec = utc_timestamp.microsecond*1000
                utc = "{} hours, {} minutes, {} seconds, {} nanoseconds".format(utc_hrs, utc_min, utc_sec, utc_msec)
                # print(utc)

                gnss_lat = float(gnss_str[2])
                lat_conv = latlon_conv(gnss_lat) 
                if gnss_str[3]=='S':
                    lat_conv = lat_conv * (-1)

                gnss_lon = float(gnss_str[4])
                lon_conv = latlon_conv(gnss_lon) 
                if gnss_str[5]=='W':
                    lon_conv = lon_conv * (-1)

                gnss_fix_quality = int(gnss_str[6])

                gnss_hdop = float(gnss_str[8])

                gnss_altitude = float(gnss_str[9])

                new_latlon = utm.from_latlon(lat_conv,lon_conv)

                # print(f'UTM_easting, UTM_northing, Zone, Letter:{new_latlon}')       
                # print('\n')
                
                msgs = gnss_msg()
                msgs.Header.stamp = rospy.Time.now()
                msgs.Header.stamp.secs = int(utc_sec)
                msgs.Header.stamp.nsecs = int(utc_msec)
                msgs.Header.frame_id = 'GNSS1_Frame'
                msgs.Latitude = lat_conv
                msgs.Longitude = lon_conv
                msgs.Altitude = gnss_altitude
                msgs.HDOP = gnss_hdop
                msgs.UTM_easting = new_latlon[0]
                msgs.UTM_northing = new_latlon[1]
                msgs.UTC = utc
                msgs.Zone = new_latlon[2]
                msgs.Letter = new_latlon[3]
                msgs.Fix_quality = gnss_fix_quality

                rospy.loginfo(msgs)
                print('\n')

                gnss_data.publish(msgs)
                rate.sleep()

if __name__==  '__main__':
    try:
        gnss_driver()

    except rospy.ROSInterruptException:
        pass                   