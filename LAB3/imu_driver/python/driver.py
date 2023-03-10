#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import serial 
import rospy
from imu_driver.msg import *
from imu_driver.msg import Vectornav
from imu_driver.srv import convert_to_quaternion, convert_to_quaternionResponse

# def e_2_q(yaw, pitch, roll):
#     ori_x = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
#     ori_y = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
#     ori_z = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
#     ori_w = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
#     return [ori_x, ori_y, ori_z, ori_w]
        
def imu_driver():
    rospy.init_node('imu')
    imu_data = rospy.Publisher('/imu', Vectornav , queue_size=15)
    rate = rospy.Rate(40) #40Hz
    seqn = 0

    connected_port = sys.argv[1]
    serial_port = rospy.get_param('~port',connected_port)
    serial_baud = rospy.get_param('~baudrate',115200)
    sampling_rate = rospy.get_param('~sampling_rate',5.0)

    sport = serial.Serial(serial_port,serial_baud,timeout = 3)
    print(serial_port)
    print('\n')
    sport.write(b'$VNWRG,07,40*XX\r')

    while not rospy.is_shutdown():
        full_imu_str = sport.readline()
        

        if full_imu_str == '':
            rospy.logwarn("No data from Vectornav IMU.")
            print("No data from Vectornav IMU")
        else:
            # if b"$VNYMR" in full_imu_str:
            #     full_imu_str = full_imu_str.decode()
            #     final=full_imu_str.split(': "')
            #     imu_str = final[1].split(',')
            #     print(imu_str)
            #     print('\n')

            if full_imu_str.startswith(b"$VNYMR") or full_imu_str.startswith(b"\r$VNYMR"):
                full_imu_str = full_imu_str.decode("utf-8")
                imu_str = full_imu_str.split(",")
                # print(imu_str)

                yaw, pitch, roll = float(imu_str[1]), float(imu_str[2]), float(imu_str[3]) 
                mag_x, mag_y, mag_z = float(imu_str[4]), float(imu_str[5]), float(imu_str[6])
                acc_x, acc_y, acc_z = float(imu_str[7]), float(imu_str[8]), float(imu_str[9])
                ang_x, ang_y, ang_z_temp = float(imu_str[10]), float(imu_str[11]), imu_str[12]
                ang_z_temps = ang_z_temp.split('*')
                ang_z = float(ang_z_temps[0])

                client = rospy.ServiceProxy("convert_to_quaternion", convert_to_quaternion)
                response = client(yaw,pitch,roll)
                (ori_x, ori_y, ori_z, ori_w) = (response.x,response.y,response.z,response.w)

                # (ori_x, ori_y, ori_z, ori_w) =e_2_q(float(np.radians(yaw)), float(np.radians(pitch)), float(np.radians(roll)))

                now = rospy.get_time()
                now_sec = int(now)
                now_nsec = int((now-now_sec)*10000000)
                seqn += 1

                msgs = Vectornav()
                msgs.raw_data = str(imu_str)

                ## Header Data Published
                msgs.Header.seq = int(seqn)
                msgs.Header.stamp = rospy.Time.now()
                msgs.Header.stamp.secs = int(now_sec)
                msgs.Header.stamp.nsecs = int(now_nsec)
                msgs.Header.frame_id = "imu1_frame"

                ## Orientation Data Published
                msgs.imu.orientation.x = ori_x
                msgs.imu.orientation.y = ori_y
                msgs.imu.orientation.z = ori_z
                msgs.imu.orientation.w = ori_w
                
                ## Acceleration Data Published
                msgs.imu.linear_acceleration.x = acc_x
                msgs.imu.linear_acceleration.y = acc_y
                msgs.imu.linear_acceleration.z = acc_z

                ## Angular Velocity Data Published
                msgs.imu.angular_velocity.x = ang_x
                msgs.imu.angular_velocity.y = ang_y
                msgs.imu.angular_velocity.z = ang_z
                
                ## Magnetic Field Data Published
                msgs.mag_field.magnetic_field.x = mag_x
                msgs.mag_field.magnetic_field.y = mag_y
                msgs.mag_field.magnetic_field.z = mag_z

                imu_data.publish(msgs)
                # rospy.loginfo(msgs)
                # print('\n')
                rate.sleep()

if __name__==  '__main__':
    try:
        imu_driver()

    except rospy.ROSInterruptException:
        pass         

    except serial.serialutil.SerialException:
        rospy.loginfo("Shutting down node...")        