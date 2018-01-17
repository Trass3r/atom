# Copyright 2015 Adobe Systems Incorporated
# All Rights Reserved.

path = require 'path'

vcInstallDir = "C:\\Program Files (x86)\\Microsoft Visual Studio 14.0\\VC"
kitsDir = "C:\\Program Files (x86)\\Windows Kits\\10"
win10ver = "10.0.16299.0"

INCLUDE = [
  (path.join vcInstallDir, "INCLUDE"),
  (path.join vcInstallDir, "ATLMFC", "INCLUDE"),
  (path.join kitsDir, "include", win10ver, "ucrt"),
  (path.join kitsDir, "include", win10ver, "shared"),
  (path.join kitsDir, "include", win10ver, "um")
]

LIBPATH = [
  (path.join vcInstallDir, "lib"),
  (path.join vcInstallDir, "atlmfc", "lib")
]

LIB = [
  (path.join vcInstallDir, "lib", "amd64"),
  (path.join vcInstallDir, "atlmfc", "lib", "amd64")
  (path.join kitsDir, "lib", win10ver, "ucrt", "x64")
  (path.join kitsDir, "lib", win10ver, "um", "x64")
]

binPath = path.join vcInstallDir, "bin", "amd64"

Path = [
  binPath,
  process.env.Path
]

module.exports.makeEnv = () ->
  INCLUDE: INCLUDE.join ";"
  LIBPATH: LIBPATH.join ";"
  LIB: LIB.join ";"
  Path: Path.join ";"

module.exports.toolPath = (basename) ->
  path.join vcInstallDir, "BIN", basename
