# C4-Real-time-pedestrian-detection

Real-Time Human Detection Using Contour Cues

based on Jianxin Wu's work

https://sites.google.com/site/wujx2001/home/c4


Real-Time Human Detection Using Contour Cues [pdf]
Jianxin Wu, Christopher Geyer, and James M. Rehg
In: Proc. The 2011 IEEE Int'l Conference on Robotics and Automation (ICRA 2011), Shanghai, China, May 2011, pp. 860-867.



**Additions on original code**

- SVM model file converted from text to binary by this change loading time is seriously reduced.

- Test interface improved.

**usage information:**

     usage:
     TestC4Detector <video_file>
     
     (note: files combined.txt.model and combined.txt.model_ must be at the same directory with executable)
     
     keys:
     space : toggle using simple post-process (NMS, non-maximal suppression)
     0     : waits to process next frame until a key pressed
     1     : doesn't wait to process next frame
     2     : resize frames 1/2
     3     : don't resize frames
     4     : resize frames 1/4

**Sample output**

<img src = "https://raw.githubusercontent.com/sturkmen72/C4-Real-time-pedestrian-detection/master/sample_output.jpg"/>

**Development objectives**

- Increased performance by using for loop in parallel.
