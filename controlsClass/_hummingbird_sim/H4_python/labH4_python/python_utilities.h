void updateKmValue() {
  if (Serial.available() > 0) {
    // There's data incoming
    Serial.print("Trying to read");
    byte x = Serial.read();
    if (x == startMarker) {
      bytesRecvd = 0;
      inProgress = true;
    }

    if(inProgress) {
      tempBuffer[bytesRecvd] = x;
      bytesRecvd ++;
    }

    if (x == endMarker) {
      inProgress = false;      
      String stringRecvd = "";
      for (byte n=1; n < bytesRecvd - 1 ; n++) {
        stringRecvd += (char)tempBuffer[n];
      }
      float kmValue = stringRecvd.toFloat();
      P.km = kmValue;
    }
  }
}
