#!/usr/bin/env python3
"""
Installation script for the clonokinetix package.
Run this script to install the package in development mode.
"""

import subprocess
import sys
import os

def install_package():
    """Install the clonokinetix package in development mode."""
    
    # Change to the project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    
    print("Installing clonokinetix package in development mode...")
    
    try:
        # Install the package in development mode
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "."])
        print("✅ Package installed successfully!")
        
        # Verify installation
        print("\nVerifying installation...")
        subprocess.check_call([sys.executable, "-c", "import clonokinetix; print('✅ clonokinetix imported successfully')"])
        
        print("\n🎉 Installation complete! You can now import clonokinetix in your notebooks.")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Installation failed: {e}")
        return False
    
    return True

if __name__ == "__main__":
    install_package()
