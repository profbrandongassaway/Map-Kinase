import base64
import os

WARNING_ICON_BASE64 = """
iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADYSURBVDhP7ZJBDsIwEEX3NtwB6AAdQN2BX+C5g7oD3cFzB3UHvYCHiF1iwGZmktz/SfIkk6RNSynAFsuyyLKU3vOcJ7Isi4wxYRh+v2+2bdN8l8V5nuf5eZ7f7/t93/e1yLIsv3fYxWmaZpvOef4N4s2jRJuW8mUoF6eZ/v7+dRGLxR3r6xez2WyOj4+Pj4+P7ff7+Xz+4+Pj6enp/X6//X6/3W6XzWYL13WmaZrOeaIoiiKfz8dxHMeZpmmaZjabdV3XbreLRCIR8F8KzwAAAABJRU5ErkJggg==
"""

# Decode and save as PNG
try:
  img_data = base64.b64decode(WARNING_ICON_BASE64)
  with open("../scripts/warning_icon.png", "wb") as f:
      f.write(img_data)
  print("PNG file saved successfully. Check 'warning_icon.png'.")
except Exception as e:
  print(f"Error decoding base64: {e}")