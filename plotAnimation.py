import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()

file = open("run1.txt")

data = file.readline()

lines = data.split(";")

x = []
y = []
index = 0
speed = 1

def animatDisplay(i):
    global index
    nums = lines[index].split(",")
    while "" in nums: nums.remove("")
    nums = [float(f) for f in nums]

    x = [i for i in range(len(nums))]
    y = nums
    ax.clear()
    ax.plot(x,y)
    ax.set_xlim([0,len(nums)])
    ax.set_ylim([-2,2])
    index += speed
    if (index >= len(lines)):
        index = 0

print("num of lines : ")
print(len(lines))
print("start animation")


ani = FuncAnimation(fig, animatDisplay, frames=int(len(lines)/speed), interval=0.5, repeat=True)
plt.show()