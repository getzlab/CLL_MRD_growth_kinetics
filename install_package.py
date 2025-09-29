#!/usr/bin/env python3
"""
Installation script for the clonokinetics package.
Run this script to install the package in development mode.
"""

import subprocess
import sys
import os

def install_package():
    """Install the clonokinetics package in development mode."""
    
    # Change to the project directory
    project_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_dir)
    
    print("Installing clonokinetics package in development mode...")
    
    try:
        # Install the package in development mode
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", "."])
        print("✅ Package installed successfully!")
        
        # Verify installation
        print("\nVerifying installation...")
        subprocess.check_call([sys.executable, "-c", "import clonokinetics; print('✅ clonokinetics imported successfully')"])
        
        print("\n🎉 Installation complete! You can now import clonokinetics in your notebooks.")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Installation failed: {e}")
        return False
    
    return True

if __name__ == "__main__":
    install_package()
