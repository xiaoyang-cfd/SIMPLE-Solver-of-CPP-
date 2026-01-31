#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 用于将路径下所有的文件编码格式从GB改成UTF

"""
Batch File Encoding Conversion Script
Function: Convert all files in specified directory from GB2312 encoding to UTF-8 encoding
"""

import os
import sys
import codecs
import argparse
from pathlib import Path

def setup_console_encoding():
    """Setup console encoding for proper output"""
    if sys.platform == 'win32':
        import ctypes
        # Windows system set console to UTF-8 encoding
        try:
            kernel32 = ctypes.windll.kernel32
            kernel32.SetConsoleOutputCP(65001)  # UTF-8
            kernel32.SetConsoleCP(65001)  # UTF-8
        except:
            pass

def convert_file_encoding(file_path, source_encoding='gb2312', target_encoding='utf-8'):
    """
    Convert encoding of a single file
    
    Args:
        file_path: Path to the file
        source_encoding: Source file encoding (default: gb2312)
        target_encoding: Target encoding (default: utf-8)
    
    Returns:
        bool: Whether conversion was successful
    """
    try:
        # First backup original file content
        try:
            with open(file_path, 'rb') as f:
                original_content = f.read()
        except Exception as e:
            print(f"  Error: Cannot read file '{file_path}': {e}")
            return False
        
        # Try to read file content with source encoding
        try:
            with codecs.open(file_path, 'r', encoding=source_encoding, errors='strict') as f:
                content = f.read()
            print(f"  Successfully read with {source_encoding} encoding")
            
        except UnicodeDecodeError:
            # If GB2312 fails, try other common encodings
            encodings_to_try = ['gbk', 'gb18030', 'cp936', 'latin-1', 'iso-8859-1', 'utf-8']
            for enc in encodings_to_try:
                try:
                    content = original_content.decode(enc)
                    print(f"  Detected encoding might be {enc}, using this encoding")
                    source_encoding = enc
                    break
                except UnicodeDecodeError:
                    continue
            else:
                print(f"  Warning: Cannot detect file encoding, using error replacement mode")
                content = original_content.decode(source_encoding, errors='replace')
        
        # Write file content with target encoding
        try:
            with codecs.open(file_path, 'w', encoding=target_encoding, errors='strict') as f:
                f.write(content)
            print(f"  Successfully converted to {target_encoding}")
            return True
            
        except Exception as e:
            print(f"  Error: Failed to write file: {e}")
            # Restore original content
            with open(file_path, 'wb') as f:
                f.write(original_content)
            return False
    
    except Exception as e:
        print(f"  Error: Unknown error processing file '{file_path}': {e}")
        return False

def batch_convert_encoding(directory, source_encoding='gb2312', target_encoding='utf-8', 
                          file_extensions=None, recursive=False):
    """
    Batch convert file encoding in directory
    
    Args:
        directory: Target directory path
        source_encoding: Source encoding (default: gb2312)
        target_encoding: Target encoding (default: utf-8)
        file_extensions: List of file extensions to process, None for all files
        recursive: Whether to process subdirectories recursively (default: False)
    """
    try:
        # Ensure directory path is string
        if isinstance(directory, Path):
            directory = str(directory)
        
        directory = os.path.abspath(directory)
        
        if not os.path.exists(directory):
            print(f"Error: Directory '{directory}' does not exist!")
            print(f"Current working directory: {os.getcwd()}")
            return
        
        print(f"Current working directory: {os.getcwd()}")
        print(f"Target directory: {directory}")
        print(f"Source encoding: {source_encoding}")
        print(f"Target encoding: {target_encoding}")
        print(f"Recursive: {'Yes' if recursive else 'No'}")
        
        if file_extensions:
            file_extensions = [ext.lower() if ext.startswith('.') else f'.{ext.lower()}'
                              for ext in file_extensions]
            print(f"File extensions filter: {file_extensions}")
        
        # Collect files to process
        files_to_process = []
        
        if recursive:
            # Process subdirectories recursively
            for root, dirs, files in os.walk(directory):
                for file in files:
                    file_path = os.path.join(root, file)
                    if not file_extensions or os.path.splitext(file)[1].lower() in file_extensions:
                        files_to_process.append(file_path)
        else:
            # Process only current directory
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isfile(item_path):
                    if not file_extensions or os.path.splitext(item)[1].lower() in file_extensions:
                        files_to_process.append(item_path)
        
        if not files_to_process:
            print(f"No matching files found in directory '{directory}'")
            return
        
        print(f"Found {len(files_to_process)} file(s) to process")
        print("-" * 60)
        
        success_count = 0
        fail_count = 0
        
        for i, file_path in enumerate(files_to_process, 1):
            rel_path = os.path.relpath(file_path, directory)
            print(f"[{i}/{len(files_to_process)}] Processing: {rel_path}")
            
            if convert_file_encoding(file_path, source_encoding, target_encoding):
                success_count += 1
            else:
                fail_count += 1
            print()  # Empty line between files
        
        print("-" * 60)
        print(f"Conversion complete!")
        print(f"Successful: {success_count} file(s)")
        print(f"Failed: {fail_count} file(s)")
        
        if fail_count > 0:
            print(f"\nNote: Some files failed to convert. They might be:")
            print("  1. Binary files (not text files)")
            print("  2. Already in UTF-8 encoding")
            print("  3. Using different encoding than specified")
            print("  4. Protected or in use by another program")
    
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user")
        return
    except Exception as e:
        print(f"\nUnexpected error: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='Batch convert file encoding (default: GB2312 -> UTF-8)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert all files in current directory
  python convert_encoding.py .
  
  # Convert files in specific directory
  python convert_encoding.py /path/to/directory
  
  # Recursively convert all subdirectories
  python convert_encoding.py . -r
  
  # Convert only specific file extensions
  python convert_encoding.py . -e txt csv ini
  
  # Specify different encoding conversion
  python convert_encoding.py . -s gbk -t utf-8
  
  # Test run without modifying files
  python convert_encoding.py . --dry-run
        """
    )
    
    parser.add_argument('directory', nargs='?', default='.',
                       help='Directory to process (default: current directory)')
    
    parser.add_argument('-s', '--source', default='gb2312',
                       help='Source file encoding (default: gb2312)')
    
    parser.add_argument('-t', '--target', default='utf-8',
                       help='Target encoding (default: utf-8)')
    
    parser.add_argument('-e', '--extensions', nargs='+',
                       help='File extensions to process (e.g., txt csv ini)')
    
    parser.add_argument('-r', '--recursive', action='store_true',
                       help='Process subdirectories recursively')
    
    parser.add_argument('--dry-run', action='store_true',
                       help='Test run without actually modifying files')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed information for each file')
    
    args = parser.parse_args()
    
    # Setup console encoding for Windows
    setup_console_encoding()
    
    if args.dry_run:
        print("DRY RUN MODE (no files will be modified)")
        print(f"Directory: {args.directory}")
        print(f"Source encoding: {args.source}")
        print(f"Target encoding: {args.target}")
        print(f"Recursive: {args.recursive}")
        print(f"File extensions: {args.extensions if args.extensions else 'All files'}")
        print("\nFiles that would be processed:")
        
        try:
            directory = Path(args.directory)
            if not directory.exists():
                print(f"Error: Directory '{directory}' does not exist!")
                return
            
            if args.recursive:
                files = list(directory.rglob('*'))
            else:
                files = list(directory.iterdir())
            
            # Filter files
            files = [f for f in files if f.is_file()]
            if args.extensions:
                extensions = [ext.lower() if ext.startswith('.') else f'.{ext.lower()}'
                             for ext in args.extensions]
                files = [f for f in files if f.suffix.lower() in extensions]
            
            for file in files:
                try:
                    rel_path = file.relative_to(directory)
                    print(f"  {rel_path}")
                except:
                    print(f"  {file}")
            
            print(f"\nTotal: {len(files)} file(s)")
            
        except Exception as e:
            print(f"Error during dry run: {e}")
        return
    
    try:
        batch_convert_encoding(
            directory=args.directory,
            source_encoding=args.source,
            target_encoding=args.target,
            file_extensions=args.extensions,
            recursive=args.recursive
        )
    except KeyboardInterrupt:
        print("\n\nOperation cancelled by user")
    except Exception as e:
        print(f"\nProgram execution error: {e}")

def convert_directory_simple(directory_path="."):
    """
    Simplified usage: call this function to process a directory
    """
    print(f"Starting conversion for directory: {directory_path}")
    batch_convert_encoding(directory_path)

if __name__ == "__main__":
    print("=" * 60)
    print("BATCH FILE ENCODING CONVERTER")
    print("=" * 60)
    print()
    
    # Method 1: Run script directly
    # python convert_encoding.py [directory_path]
    main()
    
    # Method 2: Call from Python code
    # from convert_encoding import convert_directory_simple
    # convert_directory_simple("target_directory_path")