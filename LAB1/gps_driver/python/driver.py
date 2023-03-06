#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import utm
import serial 
import rospy
from gps_driver.msg import *
from gps_driver.msg import gps_msg
from datetime import datetime


def latlon_conv(loc):
    deg = int(loc/100)
    conv = deg+(loc/100-deg)*10/6
    return conv

def gps_driver():
    rospy.init_node('gps',anonymous=True)
    gps_data = rospy.Publisher('/gps', gps_msg, queue_size=10)
    rate = rospy.Rate(10) #10Hz

    connected_port = sys.argv[1]
    serial_port = rospy.get_param('~port',connected_port)
    serial_baud = rospy.get_param('~baudrate',4800)

    sport = serial.Serial(serial_port,serial_baud,timeout = 3)
    print(serial_port)
    print('\n')

    while not rospy.is_shutdown():
        full_gps_str = sport.readline()
        #full_gps_str = full_gps_str.decode('utf-8')

        if full_gps_str == '':
            rospy.logwarn("No data from GPS Puck.")
            print("No data from puck")
        else:
            if full_gps_str.startswith(b"\r$GPGGA") or full_gps_str.startswith(b"$GPGGA"):
                full_gps_str = full_gps_str.decode()
                gps_str = full_gps_str.split(",")
                print(gps_str)

                utc_str = gps_str[1]
                utc_timestamp = datetime.strptime(utc_str,"%H%M%S.%f")
                utc_hrs = utc_timestamp.hour
                utc_min = utc_timestamp.minute
                utc_sec = utc_timestamp.second
                utc_msec = utc_timestamp.microsecond*1000
                utc = "{} hours, {} minutes, {} seconds, {} nanoseconds".format(utc_hrs, utc_min, utc_sec, utc_msec)
                print(utc)

                gps_lat = float(gps_str[2])
                lat_conv = latlon_conv(gps_lat) 
                if gps_str[3]=='S':
                    lat_conv = lat_conv * (-1)

                gps_lon = float(gps_str[4])
                lon_conv = latlon_conv(gps_lon) 
                if gps_str[5]=='W':
                    lon_conv = lon_conv * (-1)

                gps_hdop = float(gps_str[8])

                gps_altitude = float(gps_str[9])

                new_latlon = utm.from_latlon(lat_conv,lon_conv)

                print(f'UTM_easting, UTM_northing, Zone, Letter:{new_latlon}')       
                print('\n')
                
                msgs = gps_msg()

                now = rospy.get_time()
                msgs.Header.Header = "GPS_Driver"
                msgs.Header.stamp = rospy.Time.now
                msgs.Header.stamp.secs = int(utc_sec)
                msgs.Header.stamp.nsecs = int(utc_msec)
                msgs.Header.frame_id = 'GPS1_Frame'
                msgs.Latitude = lat_conv
                msgs.Longitude = lon_conv
                msgs.Altitude = gps_altitude
                msgs.HDOP = gps_hdop
                msgs.UTM_easting = new_latlon[0]
                msgs.UTM_northing = new_latlon[1]
                msgs.UTC = utc
                msgs.Zone = new_latlon[2]
                msgs.Letter = new_latlon[3]

                rospy.loginfo(msgs)
                print('\n')

                gps_data.publish(msgs)
                rate.sleep()


if __name__==  '__main__':
    try:
        gps_driver()

    except rospy.ROSInterruptException:
        pass                   