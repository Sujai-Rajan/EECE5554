#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import rospy
# from geometry_msgs.msg import Quaternion
from imu_driver.srv import convert_to_quaternion, convert_to_quaternionResponse

def euler_to_quaternion(yaw, pitch, roll):

        qx = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
        qy = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
        qz = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)
        qw = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
        # print("e2q service used")
        return [qx, qy, qz, qw]

def handle_convert_to_quaternion(req):
    
    roll = req.roll
    yaw = req.yaw
    pitch = req.pitch
    x,y,z,w = euler_to_quaternion(yaw,pitch,roll)
    # print(x,y,z,w)
    return convert_to_quaternionResponse(x,y,z,w)

def convert_to_quaternion_server():
    rospy.init_node('convert_to_quaternion_service')
    service = rospy.Service('convert_to_quaternion', convert_to_quaternion, handle_convert_to_quaternion)
    print("Ready to convert Euler angles to quaternions...\n")
    rospy.spin()

if __name__ == "__main__":
    convert_to_quaternion_server()
