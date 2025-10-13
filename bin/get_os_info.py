#!/usr/bin/env python3
import os

def get_system_id():
    info = {}
    with open("/etc/os-release", "r") as f:
        for line in f:
            if "=" in line:
                key, value = line.strip().split("=", 1)
                info[key] = value.strip('"')
    return info.get("ID", "").lower()  # 'centos', 'ubuntu'

system_id = get_system_id()

if system_id == "centos":
    cmdir = "x86_64vc"
elif system_id == "ubuntu":
    cmdir = "x86_64"
else:
    cmdir = "unknown"

print(f"Detected system: {system_id}")
print(f"cmdir = {cmdir}")

