import serial
import tkinter as tk
from threading import Thread

# Set up the serial connection (Adjust the COM port as needed)


#FIX THIS - Make sure you change 'COM8' to the correct COM port for your computer
ser = serial.Serial('/dev/cu.usbserial-AB0L9IY6', 9600)# , timeout=1)

startMarker = 125
endMarker = 126
specialByte = 124

# ----------------Km Slider Functions--------------------

def decrease_km():
    current_value = slider.get()
    new_value = max(current_value - 0.001, 0.2)
    slider.set(new_value)
    km_val.set(str(new_value))

def slider_update_km(value):
    if value == "0.000":
        value = "0.001"
    slider.set(value)
    km_val.set(value)
    km_entry.delete(0,tk.END)
    km_entry.insert(0,value)

def increase_km():
    current_value = slider.get()
    new_value = min(current_value + 0.001, 1)
    slider.set(new_value)
    km_val.set(str(new_value))

# ----------------Km Button Functions--------------------

def km_update():
    try:
        value = float(km_entry.get())
        clamped_value = max(min(value, 1), 0.001)
        km_val.set(str(clamped_value))
        slider.set(clamped_value)
    except ValueError:
        km_entry.delete(0, tk.END)
        km_entry.insert(0, 'Invalid input')


# ----------------Serial Communication Functions--------------------
            
def send_km():
    global startMarker, endMarker
    km_value = km_entry.get()
    km_adjSend = chr(startMarker) + km_value + chr(endMarker)
    byte_km = km_adjSend.encode()
    ser.write(byte_km)

def schedule_send_km():
    send_km()
    # Schedule this function to be called again after 1 milliseconds (1 second)
    root.after(100, schedule_send_km)

def read_from_arduino():
    while True:
        if ser.inWaiting() > 0:
            # serial_line = ser.readline().decode('utf-8').rstrip()  # Read the newline terminated string
            serial_line = ser.readline()
            print()          
            print(serial_line)

# ------------------------------------------------------------------------------------------------------------
            
# ----------------Start Main--------------------

# Start the thread that will read from Arduino
thread = Thread(target=read_from_arduino)
thread.daemon = True  # This thread dies when the main thread (program) exits
thread.start()

# ----------------Create the GUI--------------------

# Create the main window
root = tk.Tk()
root.title("KM Value Updater")

# Use a StringVar to observe changes
km_val = tk.StringVar(root)
km_val.set('0.338')

# ----------------Km Slider Frame--------------------

# A frame for the slider
km_slider_frame = tk.Frame(root, borderwidth=25)
km_slider_frame.pack()

# Create a "Minus" button to decrease the slider value
minus_button = tk.Button(km_slider_frame, text="-", command=decrease_km)
minus_button.pack(side=tk.LEFT)

# Create a slider
slider = tk.Scale(km_slider_frame, length=200, from_=0.1, to=0.5, resolution=0.001, orient='horizontal', command=slider_update_km)
slider.set(0.338)
slider.pack(side=tk.LEFT)

# Create a "Plus" button to increase the slider value
plus_button = tk.Button(km_slider_frame, text="+", command=increase_km)
plus_button.pack(side=tk.LEFT)

# ----------------Km Button Frame--------------------

# A frame for the button
km_button_frame = tk.Frame(root, borderwidth=25)
km_button_frame.pack()

# Create a new entry widget for inputing km
km_entry = tk.Entry(km_button_frame, width=10)
km_entry.pack(side=tk.LEFT)

# Create a button for updating Km
km_button = tk.Button(km_button_frame, text="Update Km", command=km_update)
km_button.pack(side=tk.RIGHT)

# ----------------Theta Values--------------------


# Start the periodic sending of km values
schedule_send_km()

# Run the main loop
root.mainloop()

# Close the serial connection when the GUI is closed
ser.close()
