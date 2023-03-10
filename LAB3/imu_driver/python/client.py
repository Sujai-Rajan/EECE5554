import rospy
from imu_driver.srv import convert_to_quaternion, convert_to_quaternionResponse 

def e_2_q(yaw,pitch,roll):
    rospy.init_node("convert_to_quaternion_client")
    rospy.wait_for_service("convert_to_quaternion")
    rate = rospy.Rate(1)
    while not rospy.is_shutdown():
        try:
            client = rospy.ServiceProxy("convert_to_quaternion", convert_to_quaternion)
            response = client(yaw,pitch,roll)
            rospy.loginfo(response.x)
            rate.sleep()
        
        except rospy.ServiceException as e:
            print("Service call failed",e)

if __name__ =="__main__":
    e_2_q(90,45,30) 