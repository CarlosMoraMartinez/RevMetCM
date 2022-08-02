
import sys
print("Hello from python")

print(sys.argv[1])
for i in open(sys.argv[1], 'r'):
    print(i)
for i in open(sys.argv[2], 'r'):
    print(i)
        
print("Bye from python")